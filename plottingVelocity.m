data = csvread('joint_space_velocity_segment.csv');
time = data(:,1);

joint_velocity_1 = data(:,2);
figure(1);
plot(time, joint_velocity_1, '-o');
title('Velocity of Joint 1, 2 & 4');
xlabel('Time (s)');
ylabel('Velocity (degree/s)');
hold on;

joint_velocity_2 = data(:,3);
plot(time, joint_velocity_2, '-ro');

joint_velocity_4 = data(:,5);
plot(time, joint_velocity_4, '-go');
legend("velocity of joint 1", "velocity of joint 2", "velocity of joint 4");
hold off;

joint_velocity_3 = data(:,4);
figure(2);
plot(time, joint_velocity_3, '-o');
title('Velocity of Joint 3');
xlabel('Time (s)');
ylabel('Velocity (mm)');