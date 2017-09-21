#!/usr/bin/env python

# Copyright (C) 2017 Electric Movement Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *


def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1

    alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')
    a0, a1, a2, a3, a4, a5, a6 = symbols('d0:7')
    d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
    q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8')

    dh_parameters = {
        alpha0: 0,          a0: 0,      d1: 0.75,   q1: q1,
        alpha1: -pi / 2,    a1: 0.35,   d2: 0,      q2: -pi/2 + q2,
        alpha2: 0,          a2: 1.25,   d3: 0,      q3: q3,
        alpha3: -pi / 2,    a3: -0.054, d4: 1.5,    q4: q4,
        alpha4: pi / 2,     a4: 0,      d5: 0,      q5: q5,
        alpha5: -pi / 2,    a5: 0,      d6: 0,      q6: q6,
        alpha6: 0,          a6: 0,      d7: 0.303,  q7: 0,
    }

    def dh_to_homogeneous_transform(alpha, a, d, q):
        return Matrix([
            [cos(q), -sin(q), 0, a],
            [sin(q) * cos(alpha), cos(q) * cos(alpha), -sin(alpha), -sin(alpha) * d],
            [sin(q) * sin(alpha), cos(q) * sin(alpha), cos(alpha), cos(alpha) * d],
            [0, 0, 0, 1],
        ])

    T0_1 = dh_to_homogeneous_transform(alpha0, a0, d1, q1).subs(dh_parameters)
    T1_2 = dh_to_homogeneous_transform(alpha1, a1, d2, q2).subs(dh_parameters)
    T2_3 = dh_to_homogeneous_transform(alpha2, a2, d3, q3).subs(dh_parameters)
    T3_4 = dh_to_homogeneous_transform(alpha3, a3, d4, q4).subs(dh_parameters)
    T4_5 = dh_to_homogeneous_transform(alpha4, a4, d5, q5).subs(dh_parameters)
    T5_6 = dh_to_homogeneous_transform(alpha5, a5, d6, q6).subs(dh_parameters)
    T6_EE = dh_to_homogeneous_transform(alpha6, a6, d7, q7).subs(dh_parameters)

    T0_3 = T0_1 * T1_2 * T2_3
    T0_EE = T0_3 * T3_4 * T4_5 * T5_6 * T6_EE

    # Initialize service response
    joint_trajectory_list = []
    print('# poses: {}'.format(len(req.poses)))
    for x in xrange(len(req.poses)):
        joint_trajectory_point = JointTrajectoryPoint()

        pose = req.poses[x]
        orientation = pose.orientation
        (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
            [orientation.x, orientation.y, orientation.z, orientation.w]
        )

        r, p, y = symbols('r p y')
        R_x = Matrix([
            [1, 0, 0],
            [0, cos(r), -sin(r)],
            [0, sin(r), cos(r)],
        ])
        R_y = Matrix([
            [cos(p), 0, sin(p)],
            [0, 1, 0],
            [-sin(p), 0, cos(p)],
        ])
        R_z = Matrix([
            [cos(y), -sin(y), 0],
            [sin(y), cos(y), 0],
            [0, 0, 1],
        ])
        R_correction = R_z.subs(y, radians(180)) * R_y.subs(p, radians(-90))
        R_EE = (R_z * R_y * R_x * R_correction).subs({'r': roll, 'p': pitch, 'y': yaw})

        position = pose.position
        EE = Matrix([
            [position.x],
            [position.y],
            [position.z],
        ])

        WC = EE - 0.303 * R_EE[:, 2]  # 0.303 in Z direction

        # Geometric IK method
        # theta1 defines the Z rotation of WC
        theta1 = atan2(WC[1], WC[0])

        # Link length from 2 to 3, 3 to WC and 2 to WC
        def sq(n):
            return n * n

        def dist(x, y):
            return sqrt(sq(x) + sq(y))

        link_length = [1.501, dist(dist(WC[0], WC[1]) - 0.35, WC[2] - 0.75), 1.25]

        # Cosine method to find theta2 and theta3
        from functools import reduce
        import operator

        def product(i):
            return reduce(operator.mul, i, 1)

        angles = []
        for i in range(len(link_length)):
            n = sum(i != j and sq(length) or -sq(length) for j, length in enumerate(link_length))
            d = product(length for j, length in enumerate(link_length) if i != j)
            angles.append(acos(n / (2 * d)))

        theta2 = pi / 2 - angles[0] - atan2(WC[2] - 0.75, dist(WC[0], WC[1]) - 0.35)
        theta3 = pi / 2 - (angles[1] + 0.036)

        # Calculate rotation from base frame to WC
        R0_3 = T0_3[0:3, 0:3].evalf(subs={q1: theta1, q2: theta2, q3: theta3})
        R0_3 = T0_1[0:3, 0:3] * T1_2[0:3, 0:3] * T2_3[0:3, 0:3]
        R0_3 = R0_3.evalf(subs={q1: theta1, q2: theta2, q3: theta3})

        # Calculate R3_6
        R3_6 = R0_3.inv('LU') * R_EE

        # Calculate Euler angles
        theta4 = atan2(R3_6[2, 2], -R3_6[0, 2])
        theta5 = atan2(dist(R3_6[0, 2], R3_6[2, 2]), R3_6[1, 2])
        theta6 = atan2(-R3_6[1, 1], R3_6[1, 0])

        # FK for error calculation
        FK = T0_EE.evalf(subs={q1: theta1, q2: theta2, q3: theta3, q4: theta4, q5: theta5, q6: theta6})
        your_wc = [WC[0], WC[1], WC[2]]
        your_ee = [FK[0, 3], FK[1, 3], FK[2, 3]]

        # Populate response for the IK request
        joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
        print(joint_trajectory_point)
        joint_trajectory_list.append(joint_trajectory_point)

    rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
    return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()
