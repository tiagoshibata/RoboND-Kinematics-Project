#!/usr/bin/env python

# Copyright (C) 2017 Electric Movement Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *

CALCULATE_IK_ERROR = False


def sq(n):
    return n * n


def dist(*coords):
    return sqrt(sum(sq(x) for x in coords))


def angle_between(a, b):
    angle = (a - b) % (2 * pi)
    if angle > pi:
        angle -= 2 * pi
    return angle


class IK:
    def __init__(self):
        self._initialize_direct_kinematics()
        self._initialize_inverse_kinematics()

    def _initialize_inverse_kinematics(self):
        r, p, y = symbols('r p y')

        # Rotation matrices
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
        self.R_EE = lambdify((r, p, y), R_z * R_y * R_x * R_correction)

    def _initialize_direct_kinematics(self):
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
        self.R0_3 = lambdify((q1, q2, q3), T0_3[0:3, 0:3])
        self.T0_EE = lambdify((q1, q2, q3, q4, q5, q6), (T0_3 * T3_4 * T4_5 * T5_6 * T6_EE))

    def cosine_method_angles(self, sides):
        assert len(sides) == 3
        from functools import reduce
        import operator

        def product(i):
            return reduce(operator.mul, i, 1)

        angles = []
        for i in range(len(sides)):
            n = sum(i != j and sq(length) or -sq(length) for j, length in enumerate(sides))
            d = product(length for j, length in enumerate(sides) if i != j)
            angles.append(acos(n / (2 * d)))
        return angles

    def inverse_kinematics(self, pose, last_pose=None):
        orientation = pose.orientation
        (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
            [orientation.x, orientation.y, orientation.z, orientation.w]
        )

        R_EE = Matrix(self.R_EE(roll, pitch, yaw))

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
        xy_joint_2_WC_distance = dist(WC[0], WC[1]) - 0.35
        link_length = [1.501, dist(xy_joint_2_WC_distance, WC[2] - 0.75), 1.25]

        # Cosine method to find theta2 and theta3
        angles = self.cosine_method_angles(link_length)
        # From the multiple possible solutions, the one with joint 3 above joint
        # 2 was chosen due to physical constrains of the robot
        theta2 = pi / 2 - angles[0] - atan2(WC[2] - 0.75, xy_joint_2_WC_distance)
        theta3 = pi / 2 - (angles[1] + 0.036)

        # Calculate rotation from base frame to WC
        R0_3 = self.R0_3(theta1, theta2, theta3)

        # Calculate R3_6 through inverse of R0_3 times the tranform to the end effector
        R3_6 = R0_3.transpose() * R_EE

        # Calculate Euler angles
        # Symbolic R3_6:
        # [
        # [-sin(q4) * sin(q6) + cos(q4) * cos(q5) * cos(q6),    -sin(q4) * sin(q6) + cos(q4) * cos(q5) * sin(q6),    -cos(q4) * sin(q5)],
        # [sin(q5) * cos(q6),                                   -sin(q5) * sin(q6),                                  cos(q5)],
        # [-cos(q4) * sin(q6) - sin(q4) * cos(q5) * cos(q6),    -cos(q4) * cos(q6) + sin(q4) * cos(q5) * sin(q6),    sin(q4) * sin(q5)]
        # ]
        # q4 can be derived from atan2(R3_6[2, 2], -R3_6[0, 2]) = atan2(sin(q4) * sin(q5), -(-cos(q4) * sin(q5))
        # q5 can be derived from atan2(dist(R3_6[0, 2], R3_6[2, 2]), R3_6[1, 2]) = atan2(dist(-cos(q4) * sin(q5), sin(q4) * sin(q5)), cos(q5)) = atan2(sqrt(sin(q5) ** 2 * (sin(q4) ** 2 + cos(q4) ** 2)), cos(q5)) = atan2(sin(q5), cos(q5))
        # q6 can be derived from atan2(-R3_6[1, 1], R3_6[1, 0]) = atan2(-(-sin(q5) * sin(q6)), sin(q5) * cos(q6))
        theta4 = atan2(R3_6[2, 2], -R3_6[0, 2])
        theta5 = atan2(dist(R3_6[0, 2], R3_6[2, 2]), R3_6[1, 2])
        theta6 = atan2(-R3_6[1, 1], R3_6[1, 0])
        if last_pose:
            def minimize_rotation(previous, desired):
                return min(desired, desired + pi, key=lambda x: abs(angle_between(desired, previous)))

            # If last position is known, minimize EE rotation using the fact
            # it's symmetric
            theta6 = minimize_rotation(last_pose[5], theta6)
            optimal_theta4 = minimize_rotation(last_pose[3], theta4)
            if theta4 != optimal_theta4:
                theta4 = optimal_theta4
                theta5 = -theta5

        return theta1, theta2, theta3, theta4, theta5, theta6

    def forward_kinematics(self, thetas):
        EE = self.T0_EE(*thetas)
        return EE[0:3, 3]


ik = IK()


def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print("No valid poses received")
        return -1

    # Initialize service response
    joint_trajectory_list = []
    pose = None
    for x in xrange(len(req.poses)):
        joint_trajectory_point = JointTrajectoryPoint()

        next_pose = req.poses[x]
        pose = ik.inverse_kinematics(next_pose, pose)

        # FK for error calculation
        if CALCULATE_IK_ERROR:
            calculated = ik.forward_kinematics(pose)
            error = dist(next_pose.position.x - calculated[0],
                         next_pose.position.y - calculated[1],
                         next_pose.position.z - calculated[2])
            print("End effector offset is %04.8f" % error)

        # Populate response for the IK request
        joint_trajectory_point.positions = pose
        joint_trajectory_list.append(joint_trajectory_point)

    rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
    return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print("Ready to receive an IK request")
    rospy.spin()


if __name__ == "__main__":
    IK_server()
