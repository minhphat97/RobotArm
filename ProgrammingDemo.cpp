// ProgrammingDemo.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdio.h>
#include <conio.h>
#include "ensc-488.h"
#include <iostream>
#include <string>
#include <math.h>
#include <stdlib.h> 
#include <fstream>
#include <Windows.h>
#include <algorithm>
#include <stdlib.h> 
#include <iomanip>
#include <stdlib.h>

using namespace std;
double pi = 3.14159265359;
double a2 = 195;
double a3 = 142;
double d2 = 0;
double d4 = 410;
double dT = 140;
double d0 = 475;
void printMatrix(double** Array, int size)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			cout << Array[i][j] << " ";
		}
		cout << endl;
	}
}

double** TMULT(double** brela, double** crelb)
{
	double** crela = new double*[4];
	for (int i = 0; i < 4; i++)
	{
		crela[i] = new double[4];
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			crela[i][j] = 0;
			for (int k = 0; k < 4; k++)
			{
				crela[i][j] = crela[i][j] + (brela[i][k] * crelb[k][j]);
			}
		}
	}
	return crela;
}

double** TINVERT(double** brela)
{
	double** arelb = new double*[4];
	for (int i = 0; i < 4; i++)
	{
		arelb[i] = new double[4];
	}
	arelb[0][0] = brela[0][0];
	arelb[0][1] = brela[1][0];
	arelb[0][2] = brela[2][0];
	arelb[1][0] = brela[0][1];
	arelb[1][1] = brela[1][1];
	arelb[1][2] = brela[2][1];
	arelb[2][0] = brela[0][2];
	arelb[2][1] = brela[1][2];
	arelb[2][2] = brela[2][2];
	double a = -(brela[0][0] * brela[0][3] + brela[0][1] * brela[1][3] + brela[0][2] * brela[2][3]);
	double b = -(brela[1][0] * brela[0][3] + brela[1][1] * brela[1][3] + brela[1][2] * brela[2][3]);
	double c = -(brela[2][0] * brela[0][3] + brela[2][1] * brela[1][3] + brela[2][2] * brela[2][3]);
	arelb[0][3] = a;
	arelb[1][3] = b;
	arelb[2][3] = c;
	arelb[3][0] = 0;
	arelb[3][1] = 0;
	arelb[3][2] = 0;
	arelb[3][3] = 1;
	return arelb;
}

double** KIN(double theta1, double theta2, double d3, double theta4) //theta is in degree
{
	theta1 = theta1 * (pi / 180);
	theta2 = theta2 * (pi / 180);
	theta4 = theta4 * (pi / 180);
	double** wrelb = new double* [4];
	for (int i = 0; i < 4; i++)
	{
		wrelb[i] = new double[4];
	}
	wrelb[0][0] = cos(theta1 + theta2 - theta4);
	wrelb[0][1] = sin(theta1 + theta2 - theta4);
	wrelb[0][2] = 0;
	wrelb[0][3] = a3 * cos(theta1 + theta2) + a2 * cos(theta1);
	wrelb[1][0] = sin(theta1 + theta2 - theta4);
	wrelb[1][1] = -cos(theta1 + theta2 - theta4);
	wrelb[1][2] = 0;
	wrelb[1][3] = a3 * sin(theta1 + theta2) + a2 * sin(theta1);
	wrelb[2][0] = 0;
	wrelb[2][1] = 0;
	wrelb[2][2] = -1;
	wrelb[2][3] = -d4 - d3 + d2;
	wrelb[3][0] = 0;
	wrelb[3][1] = 0;
	wrelb[3][2] = 0;
	wrelb[3][3] = 1;
	return wrelb;
}

double* WHERE(double theta1, double theta2, double d3, double theta4) //theta is in degree
{
	if (theta1 < -150 || theta1 > 150)
	{
		cout << "ERROR: The angle of JOINT 1 is out of range limit ! This JOINT configuration is invalid !" << endl;
		exit(0);
	}
	if (theta2 < -100 || theta2 > 100)
	{
		cout << "ERROR: The angle of JOINT 2 is out of range limit ! This JOINT configuration is invalid !" << endl;
		exit(0);
	}
	if (theta4 < -160 || theta4 > 160)
	{
		cout << "ERROR: The angle of JOINT 4 is out of range limit ! This JOINT configuration is invalid !" << endl;
		exit(0);
	}
	if (d3 < -200 || d3 > -100)
	{
		cout << "ERROR: The range of JOINT 3 is out of range limit ! This JOINT configuration is invalid !" << endl;
		exit(0);
	}

	double** brels = new double* [4];
	for (int i = 0; i < 4; i++)
	{
		brels[i] = new double[4];
	}
	brels[0][0] = 1;
	brels[0][1] = 0;
	brels[0][2] = 0;
	brels[0][3] = 0;
	brels[1][0] = 0;
	brels[1][1] = 1;
	brels[1][2] = 0;
	brels[1][3] = 0;
	brels[2][0] = 0;
	brels[2][1] = 0;
	brels[2][2] = 1;
	brels[2][3] = d0;
	brels[3][0] = 0;
	brels[3][1] = 0;
	brels[3][2] = 0;
	brels[3][3] = 1;

	double** wrelb = KIN(theta1, theta2, d3, theta4);

	double** trelw = new double* [4];
	for (int i = 0; i < 4; i++)
	{
		trelw[i] = new double[4];
	}
	trelw[0][0] = 1;
	trelw[0][1] = 0;
	trelw[0][2] = 0;
	trelw[0][3] = 0;
	trelw[1][0] = 0;
	trelw[1][1] = 1;
	trelw[1][2] = 0;
	trelw[1][3] = 0;
	trelw[2][0] = 0;
	trelw[2][1] = 0;
	trelw[2][2] = 1;
	trelw[2][3] = dT;
	trelw[3][0] = 0;
	trelw[3][1] = 0;
	trelw[3][2] = 0;
	trelw[3][3] = 1;

	double** wrels = TMULT(brels, wrelb);
	double** trels = TMULT(wrels, trelw);
	double* outputWHERE = new double[4];
	outputWHERE[0] = trels[0][3];
	outputWHERE[1] = trels[1][3];
	outputWHERE[2] = trels[2][3];
	double radian = atan2(trels[0][1], trels[0][0]);
	outputWHERE[3] = radian * (180 / pi);
	return outputWHERE; //return x, y, z, phi
}

double* INVKIN(double x, double y, double z, double phi) //theta is in degree
{
	if (sqrt(pow(x, 2) + pow(y, 2)) > (a2 + a3) || sqrt(pow(x, 2) + pow(y, 2)) < (a2 - a3))
	{
		cout << "ERROR: One of points is lying outside the annulus r(out) = a2 + a3 and rin = a2 - a3 !" << endl;
		exit(0);
	}
	phi = phi * (pi / 180);

	double d3 = -d4 + d2 - z + d0 - dT;
	double d3_2nd = -d4 + d2 - z + d0 - dT;
	double theta2, theta1, theta2_2nd, theta1_2nd;
	double costheta2 = (pow(x, 2) + pow(y, 2) - pow(a3, 2) - pow(a2, 2)) / (2 * a2 * a3);
	double sintheta2 = sqrt(1 - pow(costheta2, 2));
	double sintheta2_2nd = -sqrt(1 - pow(costheta2, 2));
	theta2 = atan2(sintheta2, costheta2); // function atan2 cpp will auto calculate the min value of angle
	theta2_2nd = atan2(sintheta2_2nd, costheta2);

	double denominator = pow((a3 * cos(theta2) + a2), 2) + pow((a3 * sin(theta2)), 2);
	double denominator_2nd = pow((a3 * cos(theta2_2nd) + a2), 2) + pow((a3 * sin(theta2_2nd)), 2);

	if (denominator != 0)
	{
		double numeratorCos1  = (x * (a3 * cos(theta2) + a2)) + (y * (a3 * sin(theta2)));
		double numeratorSin1 = (y * (a3 * cos(theta2) + a2)) - (x * (a3 * sin(theta2)));
		double costheta1 = numeratorCos1 / denominator;
		double sintheta1 = numeratorSin1 / denominator;
		theta1 = atan2(sintheta1, costheta1); // function atan2 cpp will auto calculate the min value of angle
	}
	else if (denominator == 0)//need compare a2 and a3 here
	{
		cout << "ERROR: This JOINT configuration is invalid !" << endl;
		exit(0);
	}

	if (denominator_2nd != 0)
	{
		double numeratorCos1_2nd = (x * (a3 * cos(theta2_2nd) + a2)) + (y * (a3 * sin(theta2_2nd)));
		double numeratorSin1_2nd = (y * (a3 * cos(theta2_2nd) + a2)) - (x * (a3 * sin(theta2_2nd)));
		double costheta1_2nd = numeratorCos1_2nd / denominator_2nd;
		double sintheta1_2nd = numeratorSin1_2nd / denominator_2nd;
		theta1_2nd = atan2(sintheta1_2nd, costheta1_2nd); // function atan2 cpp will auto calculate the min value of angle
	}
	else if (denominator_2nd == 0)//need compare a2 and a3 here
	{
		cout << "ERROR: This JOINT configuration is invalid !" << endl;
		exit(0);
	}

	theta2 = theta2 * (180 / pi);
	theta2_2nd = theta2_2nd * (180 / pi);
	theta1 = theta1 * (180 / pi);
	theta1_2nd = theta1_2nd * (180 / pi);
	phi = phi * (180 / pi);
	double theta4 = -phi + theta1 + theta2;
	double theta4_2nd = -phi + theta1_2nd + theta2_2nd;
	double* Config = new double[8];
	Config[0] = theta1;
	Config[1] = theta2;
	Config[2] = d3;
	Config[3] = theta4;
	Config[4] = theta1_2nd;
	Config[5] = theta2_2nd;
	Config[6] = d3_2nd;
	Config[7] = theta4_2nd;
	return Config;
}

double* SOLVE(double x, double y, double z, double phi) //phi is in degree
{
	double* Config = INVKIN(x, y, z, phi);
	double theta1 = Config[0];
	double theta2 = Config[1];
	double d3 = Config[2];
	double theta4 = Config[3];
	double theta1_2nd = Config[4];
	double theta2_2nd = Config[5];
	double d3_2nd = Config[6];
	double theta4_2nd = Config[7];
	double flag1 = 0;
	double flag2 = 0;
	if (theta1 < -150 || theta1 > 150)
	{
		flag1 = 1;
		cout << "ERROR: The angle of JOINT 1 of first solution is out of range limit ! This JOINT configuration is invalid !" << endl;
	}
	if (theta2 < -100 || theta2 > 100)
	{
		flag1 = 1;
		cout << "ERROR: The angle of JOINT 2 of first solution is out of range limit ! This JOINT configuration is invalid !" << endl;
	}
	if (theta4 < -160 || theta4 > 160)
	{
		flag1 = 1;
		cout << "ERROR: The angle of JOINT 4 of first solution is out of range limit ! This JOINT configuration is invalid !" << endl;
	}
	if (d3 < -200 || d3 > -100)
	{
		flag1 = 1;
		cout << "ERROR: The range of JOINT 3 of first solution is out of range limit ! This JOINT configuration is invalid !" << endl;
	}
	if (theta1_2nd < -150 || theta1_2nd > 150)
	{
		flag2 = 1;
		cout << "ERROR: The angle of JOINT 1 of second solution is out of range limit ! This JOINT configuration is invalid !" << endl;
	}
	if (theta2_2nd < -100 || theta2_2nd > 100)
	{
		flag2 = 1;
		cout << "ERROR: The angle of JOINT 2 of second solution is out of range limit ! This JOINT configuration is invalid !" << endl;
	}
	if (theta4_2nd < -160 || theta4_2nd > 160)
	{
		flag2 = 1;
		cout << "ERROR: The angle of JOINT 4 of second solution is out of range limit ! This JOINT configuration is invalid !" << endl;
	}
	if (d3_2nd < -200 || d3_2nd > -100)
	{
		flag2 = 1;
		cout << "ERROR: The range of JOINT 3 of second solution is out of range limit ! This JOINT configuration is invalid !" << endl;
	}

	double* jointConfig = new double[10];
	jointConfig[0] = theta1;
	jointConfig[1] = theta2;
	jointConfig[2] = d3;
	jointConfig[3] = theta4;
	jointConfig[4] = theta1_2nd;
	jointConfig[5] = theta2_2nd;
	jointConfig[6] = d3_2nd;
	jointConfig[7] = theta4_2nd;
	jointConfig[8] = flag1;
	jointConfig[9] = flag2;
	return jointConfig;
}

double** CUBCOEF(JOINT& initial_joint, JOINT& final_joint, double total_path_time)
{
	double joint1_initalPosition = initial_joint[0];
	double joint1_finalPosition = final_joint[0];

	double joint2_initalPosition = initial_joint[1];
	double joint2_finalPosition = final_joint[1];
	
	double joint3_initalPosition = initial_joint[2];
	double joint3_finalPosition = final_joint[2];

	double joint4_initalPosition = initial_joint[3];
	double joint4_finalPosition = final_joint[3];
	
	double a0_joint1 = joint1_initalPosition;
	double a1_joint1 = 0;
	double a2_joint1 = (3 / pow(total_path_time, 2)) * (joint1_finalPosition - joint1_initalPosition);
	double a3_joint1 = (-2 / pow(total_path_time, 3)) * (joint1_finalPosition - joint1_initalPosition);

	double a0_joint2 = joint2_initalPosition;
	double a1_joint2 = 0;
	double a2_joint2 = (3 / pow(total_path_time, 2)) * (joint2_finalPosition - joint2_initalPosition);
	double a3_joint2 = (-2 / pow(total_path_time, 3)) * (joint2_finalPosition - joint2_initalPosition);

	double a0_joint3 = joint3_initalPosition;
	double a1_joint3 = 0;
	double a2_joint3 = (3 / pow(total_path_time, 2)) * (joint3_finalPosition - joint3_initalPosition);
	double a3_joint3 = (-2 / pow(total_path_time, 3)) * (joint3_finalPosition - joint3_initalPosition);

	double a0_joint4 = joint4_initalPosition;
	double a1_joint4 = 0;
	double a2_joint4 = (3 / pow(total_path_time, 2)) * (joint4_finalPosition - joint4_initalPosition);
	double a3_joint4 = (-2 / pow(total_path_time, 3)) * (joint4_finalPosition - joint4_initalPosition);

	double** joint_coefficients = new double* [4];
	for (int i = 0; i < 4; i++)
	{
		joint_coefficients[i] = new double[4];
	}

	// Coefficients of JOINT 1 is row 1
	joint_coefficients[0][0] = a0_joint1;
	joint_coefficients[0][1] = a1_joint1;
	joint_coefficients[0][2] = a2_joint1;
	joint_coefficients[0][3] = a3_joint1;

	// Coefficients of JOINT 2 is row 2
	joint_coefficients[1][0] = a0_joint2;
	joint_coefficients[1][1] = a1_joint2;
	joint_coefficients[1][2] = a2_joint2;
	joint_coefficients[1][3] = a3_joint2;

	// Coefficients of JOINT 3 is row 3
	joint_coefficients[2][0] = a0_joint3;
	joint_coefficients[2][1] = a1_joint3;
	joint_coefficients[2][2] = a2_joint3;
	joint_coefficients[2][3] = a3_joint3;

	// Coefficients of JOINT 4 is row 4
	joint_coefficients[3][0] = a0_joint4;
	joint_coefficients[3][1] = a1_joint4;
	joint_coefficients[3][2] = a2_joint4;
	joint_coefficients[3][3] = a3_joint4;

	return joint_coefficients;
}


double* PATHPLANNER(double** joint_coefficients, double time)
{
	// Position, velocity and acceleration of JOINT 1
	double a0_joint1 = joint_coefficients[0][0];
	double a1_joint1 = joint_coefficients[0][1];
	double a2_joint1 = joint_coefficients[0][2];
	double a3_joint1 = joint_coefficients[0][3];
	double joint1Position = a0_joint1 + a1_joint1 * time + a2_joint1 * pow(time, 2) + a3_joint1 * pow(time, 3);
	double joint1Velocity = a1_joint1 + 2 * a2_joint1 * time + 3 * a3_joint1 * pow(time, 2);
	double joint1Acceleration = 2 * a2_joint1 + 6 * a3_joint1 * time;

	// Position, velocity and acceleration of JOINT 2
	double a0_joint2 = joint_coefficients[1][0];
	double a1_joint2 = joint_coefficients[1][1];
	double a2_joint2 = joint_coefficients[1][2];
	double a3_joint2 = joint_coefficients[1][3];
	double joint2Position = a0_joint2 + a1_joint2 * time + a2_joint2 * pow(time, 2) + a3_joint2 * pow(time, 3);
	double joint2Velocity = a1_joint2 + 2 * a2_joint2 * time + 3 * a3_joint2 * pow(time, 2);
	double joint2Acceleration = 2 * a2_joint2 + 6 * a3_joint2 * time;

	// Position, velocity and acceleration of JOINT 3
	double a0_joint3 = joint_coefficients[2][0];
	double a1_joint3 = joint_coefficients[2][1];
	double a2_joint3 = joint_coefficients[2][2];
	double a3_joint3 = joint_coefficients[2][3];
	double joint3Position = a0_joint3 + a1_joint3 * time + a2_joint3 * pow(time, 2) + a3_joint3 * pow(time, 3);
	double joint3Velocity = a1_joint3 + 2 * a2_joint3 * time + 3 * a3_joint3 * pow(time, 2);
	double joint3Acceleration = 2 * a2_joint3 + 6 * a3_joint3 * time;

	// Position, velocity and acceleration of JOINT 4
	double a0_joint4 = joint_coefficients[3][0];
	double a1_joint4 = joint_coefficients[3][1];
	double a2_joint4 = joint_coefficients[3][2];
	double a3_joint4 = joint_coefficients[3][3];
	double joint4Position = a0_joint4 + a1_joint4 * time + a2_joint4 * pow(time, 2) + a3_joint4 * pow(time, 3);
	double joint4Velocity = a1_joint4 + 2 * a2_joint4 * time + 3 * a3_joint4 * pow(time, 2);
	double joint4Acceleration = 2 * a2_joint4 + 6 * a3_joint4 * time;

	bool flag = 0;
	if (joint1Velocity < -150 || joint1Velocity > 150)
	{
		flag = 1;
		cout << "ERROR: The velocity of joint 1 is out of range limit !" << endl;
		//exit(0);
	}
	if (joint1Acceleration < -600 || joint1Acceleration > 600)
	{
		flag = 1;
		cout << "ERROR: The acceleration of joint 1 is out of range limit !" << endl;
		//exit(0);
	}

	if (joint2Velocity < -150 || joint2Velocity > 150)
	{
		flag = 1;
		cout << "ERROR: The velocity of joint 2 is out of range limit !" << endl;
		//exit(0);
	}
	if (joint2Acceleration < -600 || joint2Acceleration > 600)
	{
		flag = 1;
		cout << "ERROR: The acceleration of joint 2 is out of range limit !" << endl;
		//exit(0);
	}

	if (joint3Velocity < -50 || joint3Velocity > 50)
	{
		flag = 1;
		cout << "ERROR: The velocity of joint 3 is out of range limit !" << endl;
		//exit(0);
	}
	if (joint3Acceleration < -200 || joint3Acceleration > 200)
	{
		flag = 1;
		cout << "ERROR: The acceleration of joint 3 is out of range limit !" << endl;
		//exit(0);
	}

	if (joint4Velocity < -150 || joint4Velocity > 150)
	{
		flag = 1;
		cout << "ERROR: The velocity of joint 4 is out of range limit !" << endl;
		//exit(0);
	}
	if (joint4Acceleration < -600 || joint4Acceleration > 600)
	{
		flag = 1;
		cout << "ERROR: The acceleration of joint 4 is out of range limit !" << endl;
		//exit(0);
	}

	double* joint_equations = new double[14];
	joint_equations[0] = joint1Position;
	joint_equations[1] = joint1Velocity;
	joint_equations[2] = joint1Acceleration;
	joint_equations[3] = joint2Position;
	joint_equations[4] = joint2Velocity;
	joint_equations[5] = joint2Acceleration;
	joint_equations[6] = joint3Position;
	joint_equations[7] = joint3Velocity;
	joint_equations[8] = joint3Acceleration;
	joint_equations[9] = joint4Position;
	joint_equations[10] = joint4Velocity;
	joint_equations[11] = joint4Acceleration;
	joint_equations[12] = time;
	joint_equations[13] = flag;
	return joint_equations;
}
double** CartesianToJoint(JOINT& current, JOINT& Via1, JOINT& Via2, JOINT& Via3, JOINT& Goal)
{
	double** Cartesian_Space = new double* [4];
	for (int i = 0; i < 4; i++)
	{
		Cartesian_Space[i] = new double[4];
	}

	double* Via1_JointSpace = SOLVE(Via1[0], Via1[1], Via1[2], Via1[3]);
	if (Via1_JointSpace[8] == 0 && Via1_JointSpace[9] == 0)
	{
		double SlnOne_Via1 = sqrt(pow((current[0] - Via1_JointSpace[0]), 2) + pow((current[1] - Via1_JointSpace[1]), 2) + pow((current[2] - Via1_JointSpace[2]), 2) + pow((current[3] - Via1_JointSpace[3]), 2));
		double SlnTwo_Via1 = sqrt(pow((current[0] - Via1_JointSpace[4]), 2) + pow((current[1] - Via1_JointSpace[5]), 2) + pow((current[2] - Via1_JointSpace[6]), 2) + pow((current[3] - Via1_JointSpace[7]), 2));

		if (SlnOne_Via1 <= SlnTwo_Via1)
		{
			Cartesian_Space[0][0] = Via1_JointSpace[0];
			Cartesian_Space[0][1] = Via1_JointSpace[1];
			Cartesian_Space[0][2] = Via1_JointSpace[2];
			Cartesian_Space[0][3] = Via1_JointSpace[3];
		}
		else
		{
			Cartesian_Space[0][0] = Via1_JointSpace[4];
			Cartesian_Space[0][1] = Via1_JointSpace[5];
			Cartesian_Space[0][2] = Via1_JointSpace[6];
			Cartesian_Space[0][3] = Via1_JointSpace[7];
		}
	}
	else if (Via1_JointSpace[8] == 0 && Via1_JointSpace[9] == 1)
	{
		Cartesian_Space[0][0] = Via1_JointSpace[0];
		Cartesian_Space[0][1] = Via1_JointSpace[1];
		Cartesian_Space[0][2] = Via1_JointSpace[2];
		Cartesian_Space[0][3] = Via1_JointSpace[3];
	}
	else if (Via1_JointSpace[8] == 1 && Via1_JointSpace[9] == 0)
	{
		Cartesian_Space[0][0] = Via1_JointSpace[4];
		Cartesian_Space[0][1] = Via1_JointSpace[5];
		Cartesian_Space[0][2] = Via1_JointSpace[6];
		Cartesian_Space[0][3] = Via1_JointSpace[7];
	}
	else if (Via1_JointSpace[8] == 1 && Via1_JointSpace[9] == 1)
	{
		cout << "Both solutions of point VIA_1 are invalid !" << endl;
		exit(0);
	}

	double* Via2_JointSpace = SOLVE(Via2[0], Via2[1], Via2[2], Via2[3]);
	if (Via2_JointSpace[8] == 0 && Via2_JointSpace[9] == 0)
	{
		double SlnOne_Via2 = sqrt(pow((current[0] - Via2_JointSpace[0]), 2) + pow((current[1] - Via2_JointSpace[1]), 2) + pow((current[2] - Via2_JointSpace[2]), 2) + pow((current[3] - Via2_JointSpace[3]), 2));
		double SlnTwo_Via2 = sqrt(pow((current[0] - Via2_JointSpace[4]), 2) + pow((current[1] - Via2_JointSpace[5]), 2) + pow((current[2] - Via2_JointSpace[6]), 2) + pow((current[3] - Via2_JointSpace[7]), 2));

		if (SlnOne_Via2 <= SlnTwo_Via2)
		{
			Cartesian_Space[1][0] = Via2_JointSpace[0];
			Cartesian_Space[1][1] = Via2_JointSpace[1];
			Cartesian_Space[1][2] = Via2_JointSpace[2];
			Cartesian_Space[1][3] = Via2_JointSpace[3];
		}
		else
		{
			Cartesian_Space[1][0] = Via2_JointSpace[4];
			Cartesian_Space[1][1] = Via2_JointSpace[5];
			Cartesian_Space[1][2] = Via2_JointSpace[6];
			Cartesian_Space[1][3] = Via2_JointSpace[7];
		}
	}
	else if (Via2_JointSpace[8] == 0 && Via2_JointSpace[9] == 1)
	{
		Cartesian_Space[1][0] = Via2_JointSpace[0];
		Cartesian_Space[1][1] = Via2_JointSpace[1];
		Cartesian_Space[1][2] = Via2_JointSpace[2];
		Cartesian_Space[1][3] = Via2_JointSpace[3];
	}
	else if (Via2_JointSpace[8] == 1 && Via2_JointSpace[9] == 0)
	{
		Cartesian_Space[1][0] = Via2_JointSpace[4];
		Cartesian_Space[1][1] = Via2_JointSpace[5];
		Cartesian_Space[1][2] = Via2_JointSpace[6];
		Cartesian_Space[1][3] = Via2_JointSpace[7];
	}
	else if (Via2_JointSpace[8] == 1 && Via2_JointSpace[9] == 1)
	{
		cout << "Both solutions of point VIA_2 are invalid !" << endl;
		exit(0);
	}

	double* Via3_JointSpace = SOLVE(Via3[0], Via3[1], Via3[2], Via3[3]);
	if (Via3_JointSpace[8] == 0 && Via3_JointSpace[9] == 0)
	{
		double SlnOne_Via3 = sqrt(pow((current[0] - Via3_JointSpace[0]), 2) + pow((current[1] - Via3_JointSpace[1]), 2) + pow((current[2] - Via3_JointSpace[2]), 2) + pow((current[3] - Via3_JointSpace[3]), 2));
		double SlnTwo_Via3 = sqrt(pow((current[0] - Via3_JointSpace[4]), 2) + pow((current[1] - Via3_JointSpace[5]), 2) + pow((current[2] - Via3_JointSpace[6]), 2) + pow((current[3] - Via3_JointSpace[7]), 2));

		if (SlnOne_Via3 <= SlnTwo_Via3)
		{
			Cartesian_Space[2][0] = Via3_JointSpace[0];
			Cartesian_Space[2][1] = Via3_JointSpace[1];
			Cartesian_Space[2][2] = Via3_JointSpace[2];
			Cartesian_Space[2][3] = Via3_JointSpace[3];
		}
		else
		{
			Cartesian_Space[2][0] = Via3_JointSpace[4];
			Cartesian_Space[2][1] = Via3_JointSpace[5];
			Cartesian_Space[2][2] = Via3_JointSpace[6];
			Cartesian_Space[2][3] = Via3_JointSpace[7];
		}
	}
	else if (Via3_JointSpace[8] == 0 && Via3_JointSpace[9] == 1)
	{
		Cartesian_Space[2][0] = Via3_JointSpace[0];
		Cartesian_Space[2][1] = Via3_JointSpace[1];
		Cartesian_Space[2][2] = Via3_JointSpace[2];
		Cartesian_Space[2][3] = Via3_JointSpace[3];
	}
	else if (Via3_JointSpace[8] == 1 && Via3_JointSpace[9] == 0)
	{
		Cartesian_Space[2][0] = Via3_JointSpace[4];
		Cartesian_Space[2][1] = Via3_JointSpace[5];
		Cartesian_Space[2][2] = Via3_JointSpace[6];
		Cartesian_Space[2][3] = Via3_JointSpace[7];
	}
	else if (Via3_JointSpace[8] == 1 && Via3_JointSpace[9] == 1)
	{
		cout << "Both solutions of point VIA_3 are invalid !" << endl;
		exit(0);
	}

	double* Goal_JointSpace = SOLVE(Goal[0], Goal[1], Goal[2], Goal[3]);
	if (Goal_JointSpace[8] == 0 && Goal_JointSpace[9] == 0)
	{
		double SlnOne_Goal = sqrt(pow((current[0] - Goal_JointSpace[0]), 2) + pow((current[1] - Goal_JointSpace[1]), 2) + pow((current[2] - Goal_JointSpace[2]), 2) + pow((current[3] - Goal_JointSpace[3]), 2));
		double SlnTwo_Goal = sqrt(pow((current[0] - Goal_JointSpace[4]), 2) + pow((current[1] - Goal_JointSpace[5]), 2) + pow((current[2] - Goal_JointSpace[6]), 2) + pow((current[3] - Goal_JointSpace[7]), 2));

		if (SlnOne_Goal <= SlnTwo_Goal)
		{
			Cartesian_Space[3][0] = Goal_JointSpace[0];
			Cartesian_Space[3][1] = Goal_JointSpace[1];
			Cartesian_Space[3][2] = Goal_JointSpace[2];
			Cartesian_Space[3][3] = Goal_JointSpace[3];
		}
		else
		{
			Cartesian_Space[3][0] = Goal_JointSpace[4];
			Cartesian_Space[3][1] = Goal_JointSpace[5];
			Cartesian_Space[3][2] = Goal_JointSpace[6];
			Cartesian_Space[3][3] = Goal_JointSpace[7];
		}
	}
	else if (Goal_JointSpace[8] == 0 && Goal_JointSpace[9] == 1)
	{
		Cartesian_Space[3][0] = Goal_JointSpace[0];
		Cartesian_Space[3][1] = Goal_JointSpace[1];
		Cartesian_Space[3][2] = Goal_JointSpace[2];
		Cartesian_Space[3][3] = Goal_JointSpace[3];
	}
	else if (Goal_JointSpace[8] == 1 && Goal_JointSpace[9] == 0)
	{
		Cartesian_Space[3][0] = Goal_JointSpace[4];
		Cartesian_Space[3][1] = Goal_JointSpace[5];
		Cartesian_Space[3][2] = Goal_JointSpace[6];
		Cartesian_Space[3][3] = Goal_JointSpace[7];
	}
	else if (Goal_JointSpace[8] == 1 && Goal_JointSpace[9] == 1)
	{
		cout << "Both solutions of point GOAL are invalid !" << endl;
		exit(0);
	}

	return Cartesian_Space;

}

double Segment_Time(JOINT& initial_point, JOINT & final_point)
{
	double final_time_joint1, final_time_joint2, final_time_joint3, final_time_joint4;

	double joint1_time1 = (3.0 / (2.0 * 150.0)) * abs(final_point[0] - initial_point[0]);
	double joint1_time2 = sqrt((6.0 / 600.0) * abs(final_point[0] - initial_point[0]));
	final_time_joint1 = max(joint1_time1, joint1_time2);
	cout << "final_time_joint1: " << final_time_joint1 << endl;

	double joint2_time1 = (3.0 / (2.0 * 150.0)) * abs(final_point[1] - initial_point[1]);
	double joint2_time2 = sqrt((6.0 / 600.0) * abs(final_point[1] - initial_point[1]));
	final_time_joint2 = max(joint2_time1, joint2_time2);
	cout << "final_time_joint2: " << final_time_joint2 << endl;

	double joint3_time1 = (3.0 / (2.0 * 50.0)) * abs(final_point[2] - initial_point[2]);
	double joint3_time2 = sqrt((6.0 / 200.0) * abs(final_point[2] - initial_point[2]));
	final_time_joint3 = max(joint3_time1, joint3_time2);
	cout << "final_time_joint3: " << final_time_joint3 << endl;

	double joint4_time1 = (3.0 / (2.0 * 150.0)) * abs(final_point[3] - initial_point[3]);
	double joint4_time2 = sqrt((6.0 / 600.0) * abs(final_point[3] - initial_point[3]));
	final_time_joint4 = max(joint4_time1, joint4_time2);
	cout << "final_time_joint4: " << final_time_joint4 << endl;

	double path_total_time = max(max(final_time_joint1, final_time_joint2), max(final_time_joint3, final_time_joint4));
	return path_total_time;
}

int main(int argc, char* argv[])
{
	JOINT q1 = { 0, 0, 0, 0 }; // Used for ForwardKin
	//JOINT q2Inverse = { 0, 0, 0, 0 }; // Used for InvreseKin 
	printf("Keep this window in focus, and...\n");
	char ch;
	int c;
	bool grabber = 0;
	const int ESC = 27;

	printf("Press any key to continue \n");
	printf("Press ESC to exit \n");

	c = _getch();

	while (1)
	{
		if (c != ESC)
		{
			printf("Press '1' for ForwardKin, '2' for InverseKin, '3' for Tracjectory Planner\n");
			ch = _getch();
			if (ch == '1')
			{
				printf("Press '1' to enter parameter values, '2' to go back\n");
				char ch1 = _getch();
				if (ch1 == '1')
				{
					printf("Enter theta1 value:\n");
					cin >> q1[0];
					printf("Enter theta2 value:\n");
					cin >> q1[1];
					printf("Enter d3 value:\n");
					cin >> q1[2];
					printf("Enter theta4 value:\n");
					cin >> q1[3];
					double* Where = WHERE(q1[0], q1[1], q1[2], q1[3]);
					printf("The resulting configuration is: \n");
					cout << "X (mm): " << Where[0] << endl;
					cout << "Y (mm): " << Where[1] << endl;
					cout << "Z (mm): " << Where[2] << endl;
					cout << "Phi angle (degree): " << Where[3] << endl;
					MoveToConfiguration(q1);
					/*cout << "To close the grabber press1, to open press0: \n";
					cin >> grabber;
					Grasp(grabber);*/
				}
				else if (ch1 == '2') {}
				else
					break;
			}
			else if (ch == '2')
			{
				printf("Press '1' to enter parameter values, '2' to go back\n");
				double x, y, z, phi;
				char ch2 = _getch();
				if (ch2 == '1')
				{
					printf("Enter the x value: \n");
					cin >> x;
					printf("Enter the y value: \n");
					cin >> y;
					printf("Enter the z value: \n");
					cin >> z;
					printf("Enter the phi value: \n");
					cin >> phi;
					double* q2Solve = SOLVE(x, y, z, phi); //(x, y, z, phi)
					cout << "theta1_1stSolution: " << q2Solve[0] << endl;
					cout << "theta2_1stSolution: " << q2Solve[1] << endl;
					cout << "d3_1stSolution: " << q2Solve[2] << endl;
					cout << "theta3_1stSolution: " << q2Solve[3] << endl;
					cout << endl;
					cout << "theta1_2ndSolution: " << q2Solve[4] << endl;
					cout << "theta2_2ndSolution: " << q2Solve[5] << endl;
					cout << "d3_2ndSolution: " << q2Solve[6] << endl;
					cout << "theta3_2ndSolution: " << q2Solve[7] << endl;
					cout << endl;
					if (q2Solve[8] == 0 && q2Solve[9] == 0)
					{
						JOINT CurrentSate;
						GetConfiguration(CurrentSate);

						double SlnOne = sqrt(pow((CurrentSate[0] - q2Solve[0]), 2) + pow((CurrentSate[1] - q2Solve[1]), 2) + pow((CurrentSate[3] - q2Solve[3]), 2));
						double SlnTwo = sqrt(pow((CurrentSate[0] - q2Solve[4]), 2) + pow((CurrentSate[1] - q2Solve[5]), 2) + pow((CurrentSate[3] - q2Solve[7]), 2));

						if (SlnOne <= SlnTwo)
						{
							JOINT q2Inverse = { q2Solve[0], q2Solve[1], q2Solve[2], q2Solve[3] };
							cout << "The L2 norm of solution one is: " << SlnOne << endl;
							cout << "The L2 norm of solution two is: " << SlnTwo << endl;
							cout << "The first solution is chosen." << endl << endl;
							MoveToConfiguration(q2Inverse);
						}
						else
						{
							JOINT q2Inverse = { q2Solve[4], q2Solve[5], q2Solve[6], q2Solve[7] };
							cout << "The L2 norm of solution one is: " << SlnOne << endl;
							cout << "The L2 norm of solution two is: " << SlnTwo << endl;
							cout << "The second solution is chosen." << endl << endl;
							MoveToConfiguration(q2Inverse);
						}
					}
					else if (q2Solve[8] == 0 && q2Solve[9] == 1)
					{
						JOINT q2Inverse = { q2Solve[0], q2Solve[1], q2Solve[2], q2Solve[3] };
						cout << "Only first solution is valid." << endl;
						MoveToConfiguration(q2Inverse);
					}

					else if (q2Solve[8] == 1 && q2Solve[9] == 0)
					{
						JOINT q2Inverse = { q2Solve[4], q2Solve[5], q2Solve[6], q2Solve[7] };
						cout << "Only second solution is valid." << endl;
						MoveToConfiguration(q2Inverse);
					}

					else if (q2Solve[8] == 1 && q2Solve[9] == 1)
					{
						cout << "Both solutions are invalid !" << endl;
						exit(0);
					}

					/*cout << "To close the grabber press1, to open press0: \n";
					cin >> grabber;
					Grasp(grabber);*/
				}
				else if (ch2 == '2') {}
				else
					break;
			}
			else if (ch == '3')
			{
				printf("Press '1' to enter parameter values, '2' to go back\n");
				char ch3 = _getch();
				if (ch3 == '1')
				{
					double x1, y1, z1, phi1;
					double x2, y2, z2, phi2;
					double x3, y3, z3, phi3;
					double xGoal, yGoal, zGoal, phiGoal;

					cout << "Please enter the x-coordinator of VIA point 1: ";
					cin >> x1;
					cout << "Please enter the y-coordinator of VIA point 1: ";
					cin >> y1;
					cout << "Please enter the z-coordinator of VIA point 1: ";
					cin >> z1;
					cout << "Please enter the Phi of VIA point 1: ";
					cin >> phi1;

					cout << "Please enter the x-coordinator of VIA point 2: ";
					cin >> x2;
					cout << "Please enter the y-coordinator of VIA point 2: ";
					cin >> y2;
					cout << "Please enter the z-coordinator of VIA point 2: ";
					cin >> z2;
					cout << "Please enter the Phi of VIA point 2: ";
					cin >> phi2;

					cout << "Please enter the x-coordinator of VIA point 3: ";
					cin >> x3;
					cout << "Please enter the y-coordinator of VIA point 3: ";
					cin >> y3;
					cout << "Please enter the z-coordinator of VIA point 3: ";
					cin >> z3;
					cout << "Please enter the Phi of VIA point 3: ";
					cin >> phi3;

					cout << "Please enter the x-coordinator of point Goal: ";
					cin >> xGoal;
					cout << "Please enter the y-coordinator of point Goal: ";
					cin >> yGoal;
					cout << "Please enter the z-coordinator of point Goal: ";
					cin >> zGoal;
					cout << "Please enter the Phi of point Goal: ";
					cin >> phiGoal;

					JOINT Via1 = { x1, y1, z1, phi1 };
					JOINT Via2 = { x2, y2, z2, phi2 };
					JOINT Via3 = { x3, y3, z3, phi3 };
					JOINT Goal = { xGoal, yGoal, zGoal, phiGoal };
					JOINT current;
					GetConfiguration(current);
					cout << "Current: " << current[0] << " " << current[1] << " " << current[2] << " " << current[3] << endl;
					double** Joint_Space = CartesianToJoint(current, Via1, Via2, Via3, Goal);
					JOINT Via1_Cartesian = { Joint_Space[0][0], Joint_Space[0][1], Joint_Space[0][2], Joint_Space[0][3] };
					JOINT Via2_Cartesian = { Joint_Space[1][0], Joint_Space[1][1], Joint_Space[1][2], Joint_Space[1][3] };
					JOINT Via3_Cartesian = { Joint_Space[2][0], Joint_Space[2][1], Joint_Space[2][2], Joint_Space[2][3] };
					JOINT Goal_Cartesian = { Joint_Space[3][0], Joint_Space[3][1], Joint_Space[3][2], Joint_Space[3][3] };

					double total_time_segment_1 = Segment_Time(current, Via1_Cartesian);
					cout << "path_total_time_segment_1: " << total_time_segment_1 << endl;

					double total_time_segment_2 = Segment_Time(Via1_Cartesian, Via2_Cartesian);
					cout << "path_total_time_segment_2: " << total_time_segment_2 << endl;

					double total_time_segment_3 = Segment_Time(Via2_Cartesian, Via3_Cartesian);
					cout << "path_total_time_segment_3: " << total_time_segment_3 << endl;

					double total_time_segment_4 = Segment_Time(Via3_Cartesian, Goal_Cartesian);
					cout << "path_total_time_segment_4: " << total_time_segment_4 << endl;

					double total_trajectory_time;
					cout << "The minimum amount of time for trajectory from Intial point to First Via Point to avoid exceeding velocity and acceleration limit is: " << total_time_segment_1 << endl;
					cout << "The minimum amount of time for trajectory from First Via Point point to Second Via Point to avoid exceeding velocity and acceleration limit is: " << total_time_segment_2 << endl;
					cout << "The minimum amount of time for trajectory from Second Via Point point to Third Via Point to avoid exceeding velocity and acceleration limit is: " << total_time_segment_3 << endl;
					cout << "The minimum amount of time for trajectory from Third Via Point point to Final Point to avoid exceeding velocity and acceleration limit is: " << total_time_segment_4 << endl;
					cout << "Therfore, we strongly recommend to choose total trajectory time bigger than " << total_time_segment_1 << " + " << total_time_segment_2 << " + " << total_time_segment_3 << " + " << total_time_segment_4 << " = " << total_time_segment_1 + total_time_segment_2 + total_time_segment_3 + total_time_segment_4 << endl;
					cout << "Please enter the total trajectory time: ";
					cin >> total_trajectory_time;
					double remaining_time = total_trajectory_time - (total_time_segment_1 + total_time_segment_2 + total_time_segment_3 + total_time_segment_4);
					total_time_segment_1 = total_time_segment_1 + remaining_time / 4;
					total_time_segment_2 = total_time_segment_2 + remaining_time / 4;
					total_time_segment_3 = total_time_segment_3 + remaining_time / 4;
					total_time_segment_4 = total_time_segment_4 + remaining_time / 4;
					cout << "total_1: " << total_time_segment_1 << endl;
					cout << "total_2: " << total_time_segment_2 << endl;
					cout << "total_3: " << total_time_segment_3 << endl;
					cout << "total_4: " << total_time_segment_4 << endl;
					cout << "time_remaining: " << remaining_time << endl;
					double** coef_segment_1 = CUBCOEF(current, Via1_Cartesian, total_time_segment_1);
					double** coef_segment_2 = CUBCOEF(Via1_Cartesian, Via2_Cartesian, total_time_segment_2);
					double** coef_segment_3 = CUBCOEF(Via2_Cartesian, Via3_Cartesian, total_time_segment_3);
					double** coef_segment_4 = CUBCOEF(Via3_Cartesian, Goal_Cartesian, total_time_segment_4);

					//=============================================================PLOTTING=====================================================================
					double plotting_time_segment1 = 0;
					double step_time_segment_1_plotting = total_time_segment_1 / 10;
					ofstream fout_position_segment("D:/ENSC488/CopytoYourDir2023/CopytoYourDir2023/joint_space_position_segment.csv");
					ofstream fout_velocity_segment("D:/ENSC488/CopytoYourDir2023/CopytoYourDir2023/joint_space_velocity_segment.csv"); 
					ofstream fout_acceleration_segment("D:/ENSC488/CopytoYourDir2023/CopytoYourDir2023/joint_space_acceleration_segment.csv");
					while (true)
					{
						double* planner = PATHPLANNER(coef_segment_1, plotting_time_segment1);
						JOINT position = { planner[0], planner[3], planner[6], planner[9] };
						JOINT velocity = { planner[1], planner[4], planner[7], planner[10] };
						JOINT acceleration = { planner[2], planner[5], planner[8], planner[11] };
						double* plotXY = WHERE(position[0], position[1], position[2], position[3]);
						fout_position_segment << plotting_time_segment1 << "," << planner[0] << "," << planner[3] << "," << planner[6] << "," << planner[9] << "," << plotXY[0] << "," << plotXY[1] << endl;
						fout_velocity_segment << plotting_time_segment1 << "," << planner[1] << "," << planner[4] << "," << planner[7] << "," << planner[10] << endl;
						fout_acceleration_segment << plotting_time_segment1 << "," << planner[2] << "," << planner[5] << "," << planner[8] << "," << planner[11] << endl;
						cout << "DONE" << endl;
						plotting_time_segment1 = plotting_time_segment1 + step_time_segment_1_plotting;
						if (plotting_time_segment1 > total_time_segment_1)
						{
							break;
						}
					}

					double plotting_time_segment2 = 0;
					double step_time_segment_2_plotting = total_time_segment_2 / 10;
					while (true)
					{
						double* planner = PATHPLANNER(coef_segment_2, plotting_time_segment2);
						JOINT position = { planner[0], planner[3], planner[6], planner[9] };
						JOINT velocity = { planner[1], planner[4], planner[7], planner[10] };
						JOINT acceleration = { planner[2], planner[5], planner[8], planner[11] };
						double* plotXY = WHERE(position[0], position[1], position[2], position[3]);
						fout_position_segment << plotting_time_segment2 + plotting_time_segment1 << "," << planner[0] << "," << planner[3] << "," << planner[6] << "," << planner[9] << "," << plotXY[0] << "," << plotXY[1] << endl;
						fout_velocity_segment << plotting_time_segment2 + plotting_time_segment1 << "," << planner[1] << "," << planner[4] << "," << planner[7] << "," << planner[10] << endl;
						fout_acceleration_segment << plotting_time_segment2 + plotting_time_segment1 << "," << planner[2] << "," << planner[5] << "," << planner[8] << "," << planner[11] << endl;
						cout << "DONE" << endl;
						plotting_time_segment2 = plotting_time_segment2 + step_time_segment_2_plotting;
						if (plotting_time_segment2 > total_time_segment_2)
						{
							break;
						}
					}
					
					double plotting_time_segment3 = 0;
					double step_time_segment_3_plotting = total_time_segment_3 / 10;
					while (true)
					{
						double* planner = PATHPLANNER(coef_segment_3, plotting_time_segment3);
						JOINT position = { planner[0], planner[3], planner[6], planner[9] };
						JOINT velocity = { planner[1], planner[4], planner[7], planner[10] };
						JOINT acceleration = { planner[2], planner[5], planner[8], planner[11] };
						double* plotXY = WHERE(position[0], position[1], position[2], position[3]);
						fout_position_segment << plotting_time_segment3 + plotting_time_segment2 + plotting_time_segment1 << "," << planner[0] << "," << planner[3] << "," << planner[6] << "," << planner[9] << "," << plotXY[0] << "," << plotXY[1] << endl;
						fout_velocity_segment << plotting_time_segment3 + plotting_time_segment2 + plotting_time_segment1 << "," << planner[1] << "," << planner[4] << "," << planner[7] << "," << planner[10] << endl;
						fout_acceleration_segment << plotting_time_segment3 + plotting_time_segment2 + plotting_time_segment1 << "," << planner[2] << "," << planner[5] << "," << planner[8] << "," << planner[11] << endl;
						cout << "DONE" << endl;
						plotting_time_segment3 = plotting_time_segment3 + step_time_segment_3_plotting;
						if (plotting_time_segment3 > total_time_segment_3)
						{
							break;
						}
					}

					double plotting_time_segment4 = 0;
					double step_time_segment_4_plotting = total_time_segment_4 / 10;
					while (true)
					{
						double* planner = PATHPLANNER(coef_segment_4, plotting_time_segment4);
						JOINT position = { planner[0], planner[3], planner[6], planner[9] };
						JOINT velocity = { planner[1], planner[4], planner[7], planner[10] };
						JOINT acceleration = { planner[2], planner[5], planner[8], planner[11] };
						double* plotXY = WHERE(position[0], position[1], position[2], position[3]);
						fout_position_segment << plotting_time_segment4 + plotting_time_segment3 + plotting_time_segment2 + plotting_time_segment1 << "," << planner[0] << "," << planner[3] << "," << planner[6] << "," << planner[9] << "," << plotXY[0] << "," << plotXY[1] << endl;
						fout_velocity_segment << plotting_time_segment4 + plotting_time_segment3 + plotting_time_segment2 + plotting_time_segment1 << "," << planner[1] << "," << planner[4] << "," << planner[7] << "," << planner[10] << endl;
						fout_acceleration_segment << plotting_time_segment4 + plotting_time_segment3 + plotting_time_segment2 + plotting_time_segment1 << "," << planner[2] << "," << planner[5] << "," << planner[8] << "," << planner[11] << endl;
						cout << "DONE" << endl;
						plotting_time_segment4 = plotting_time_segment4 + step_time_segment_4_plotting;
						if (plotting_time_segment4 > total_time_segment_4)
						{
							break;
						}
					}
					fout_position_segment.close();
					fout_velocity_segment.close();
					fout_acceleration_segment.close();


					//==============================================MOVING ROBOT ARM=====================================================
					double timing = 0;
					double step_time_segment_1 = total_time_segment_1 / 10;
					while (true)
					{
						cout << "################################################" << endl;
						cout << "At timing = " << timing << endl;
						double* planner = PATHPLANNER(coef_segment_1, timing);
						if (planner[13] == 1)
						{
							exit(0);
						}
						JOINT position = { planner[0], planner[3], planner[6], planner[9] };
						JOINT velocity = { planner[1], planner[4], planner[7], planner[10] };
						JOINT acceleration = { planner[2], planner[5], planner[8], planner[11] };
						bool check1 = MoveWithConfVelAcc(position, velocity, acceleration);
						Sleep(step_time_segment_1 * 1000);
						timing = timing + step_time_segment_1;
						if (timing > total_time_segment_1)
						{
							JOINT velocity_zero = { 0, 0, 0, 0 };
							MoveWithConfVelAcc(position, velocity_zero, acceleration);
							break;
						}
					}

					timing = 0;
					double step_time_segment_2 = total_time_segment_2 / 10;
					while (true)
					{
						cout << "################################################" << endl;
						cout << "At timing = " << timing << endl;
						double* planner = PATHPLANNER(coef_segment_2, timing);
						if (planner[13] == 1)
						{
							exit(0);
						}
						JOINT position = { planner[0], planner[3], planner[6], planner[9] };
						JOINT velocity = { planner[1], planner[4], planner[7], planner[10] };
						JOINT acceleration = { planner[2], planner[5], planner[8], planner[11] };
						bool check2 = MoveWithConfVelAcc(position, velocity, acceleration);
						Sleep(step_time_segment_2 * 1000);
						timing = timing + step_time_segment_2;
						if (timing > total_time_segment_2)
						{
							JOINT velocity_zero = { 0, 0, 0, 0 };
							MoveWithConfVelAcc(position, velocity_zero, acceleration);
							break;
						}
					}

					timing = 0;
					double step_time_segment_3 = total_time_segment_3 / 10;
					while (true)
					{
						cout << "################################################" << endl;
						cout << "At timing = " << timing << endl;
						double* planner = PATHPLANNER(coef_segment_3, timing);
						if (planner[13] == 1)
						{
							exit(0);
						}
						JOINT position = { planner[0], planner[3], planner[6], planner[9] };
						JOINT velocity = { planner[1], planner[4], planner[7], planner[10] };
						JOINT acceleration = { planner[2], planner[5], planner[8], planner[11] };
						bool check3 = MoveWithConfVelAcc(position, velocity, acceleration);
						Sleep(step_time_segment_3 * 1000);
						timing = timing + step_time_segment_3;
						if (timing > total_time_segment_3)
						{
							JOINT velocity_zero = { 0, 0, 0, 0 };
							MoveWithConfVelAcc(position, velocity_zero, acceleration);
							break;
						}
					}

					timing = 0;
					double step_time_segment_4 = total_time_segment_4 / 10;
					while (true)
					{
						cout << "################################################" << endl;
						cout << "At timing = " << timing << endl;
						double* planner = PATHPLANNER(coef_segment_4, timing);
						if (planner[13] == 1)
						{
							exit(0);
						}
						JOINT position = { planner[0], planner[3], planner[6], planner[9] };
						JOINT velocity = { planner[1], planner[4], planner[7], planner[10] };
						JOINT acceleration = { planner[2], planner[5], planner[8], planner[11] };
						bool check4 = MoveWithConfVelAcc(position, velocity, acceleration);
						Sleep(step_time_segment_4 * 1000);
						timing = timing + step_time_segment_4;
						if (timing > total_time_segment_4)
						{
							JOINT velocity_zero = { 0, 0, 0, 0 };
							MoveWithConfVelAcc(position, velocity_zero, acceleration);
							break;
						}
					}
				}
				else if (ch3 == '2') {}
				else
					break;
			}
		}
		else
			break;
	}
	return 0;
}


