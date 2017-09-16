from sympy import *
from time import time
from mpmath import radians
import tf

'''
Format of test case is [ [[EE position],[EE orientation as quaternions]],[WC location],[joint angles]]
You can generate additional test cases by setting up your kuka project and running `$ roslaunch kuka_arm forward_kinematics.launch`
From here you can adjust the joint angles to find thetas, use the gripper to extract positions and orientation (in quaternion xyzw) and lastly use link 5
to find the position of the wrist center. These newly generated test cases can be added to the test_cases dictionary.
'''

test_cases = {1:[[[2.16135,-1.42635,1.55109],
                  [0.708611,0.186356,-0.157931,0.661967]],
                  [1.89451,-1.44302,1.69366],
                  [-0.65,0.45,-0.36,0.95,0.79,0.49]],
              2:[[[-0.56754,0.93663,3.0038],
                  [0.62073, 0.48318,0.38759,0.480629]],
                  [-0.638,0.64198,2.9988],
                  [-0.79,-0.11,-2.33,1.94,1.14,-3.68]],
              3:[[[-1.3863,0.02074,0.90986],
                  [0.01735,-0.2179,0.9025,0.371016]],
                  [-1.1669,-0.17989,0.85137],
                  [-2.99,-0.12,0.94,4.06,1.29,-4.12]],
              }


def test_code(test_case):
    ## Set up code
    ## Do not modify!
    x = 0
    class Position:
        def __init__(self,EE_pos):
            self.x = EE_pos[0]
            self.y = EE_pos[1]
            self.z = EE_pos[2]
    class Orientation:
        def __init__(self,EE_ori):
            self.x = EE_ori[0]
            self.y = EE_ori[1]
            self.z = EE_ori[2]
            self.w = EE_ori[3]

    position = Position(test_case[0][0])
    orientation = Orientation(test_case[0][1])

    class Combine:
        def __init__(self,position,orientation):
            self.position = position
            self.orientation = orientation

    comb = Combine(position,orientation)

    class Pose:
        def __init__(self,comb):
            self.poses = [comb]

    req = Pose(comb)
    start_time = time()

    # IK calculation
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

    # Error analysis
    print ("\nTotal run time to calculate joint angles from pose is %04.4f seconds" % (time()-start_time))

    # Find WC error
    if not(sum(your_wc)==3):
        wc_x_e = abs(your_wc[0]-test_case[1][0])
        wc_y_e = abs(your_wc[1]-test_case[1][1])
        wc_z_e = abs(your_wc[2]-test_case[1][2])
        wc_offset = sqrt(wc_x_e**2 + wc_y_e**2 + wc_z_e**2)
        print ("\nWrist error for x position is: %04.8f" % wc_x_e)
        print ("Wrist error for y position is: %04.8f" % wc_y_e)
        print ("Wrist error for z position is: %04.8f" % wc_z_e)
        print ("Overall wrist offset is: %04.8f units" % wc_offset)

    # Find theta errors
    t_1_e = abs(theta1-test_case[2][0])
    t_2_e = abs(theta2-test_case[2][1])
    t_3_e = abs(theta3-test_case[2][2])
    t_4_e = abs(theta4-test_case[2][3])
    t_5_e = abs(theta5-test_case[2][4])
    t_6_e = abs(theta6-test_case[2][5])
    print ("\nTheta 1 error is: %04.8f" % t_1_e)
    print ("Theta 2 error is: %04.8f" % t_2_e)
    print ("Theta 3 error is: %04.8f" % t_3_e)
    print ("Theta 4 error is: %04.8f" % t_4_e)
    print ("Theta 5 error is: %04.8f" % t_5_e)
    print ("Theta 6 error is: %04.8f" % t_6_e)
    print ("\n**These theta errors may not be a correct representation of your code, due to the fact \
           \nthat the arm can have muliple positions. It is best to add your forward kinematics to \
           \nconfirm whether your code is working or not**")
    print (" ")

    # Find FK EE error
    if not(sum(your_ee)==3):
        ee_x_e = abs(your_ee[0]-test_case[0][0][0])
        ee_y_e = abs(your_ee[1]-test_case[0][0][1])
        ee_z_e = abs(your_ee[2]-test_case[0][0][2])
        ee_offset = sqrt(ee_x_e**2 + ee_y_e**2 + ee_z_e**2)
        print ("\nEnd effector error for x position is: %04.8f" % ee_x_e)
        print ("End effector error for y position is: %04.8f" % ee_y_e)
        print ("End effector error for z position is: %04.8f" % ee_z_e)
        print ("Overall end effector offset is: %04.8f units \n" % ee_offset)


if __name__ == "__main__":
    # Change test case number for different scenarios
    for test_case in test_cases.values():
        print('=' * 40)
        test_code(test_case)
