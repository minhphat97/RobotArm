data = csvread('joint_space_acceleration_segment.csv');
time = data(:,1);

joint_acceleration_1 = data(:,2);
figure(1);
plot(time, joint_acceleration_1, '-o');
title('Acceleration of Joint 1, 2 & 4');
xlabel('Time (s)');
ylabel('Acceleration (degree/s^2)');
hold on;

joint_acceleration_2 = data(:,3);
plot(time, joint_acceleration_2, '-ro');

joint_acceleration_4 = data(:,5);
plot(time, joint_acceleration_4, '-go');
legend("acceleration of joint 1", "acceleration of joint 2", "acceleration of joint 4");
hold off;

joint_acceleration_3 = data(:,4);
figure(2);
plot(time, joint_acceleration_3, '-o');
title('Acceleration of Joint 3');
xlabel('Time (s)');
ylabel('Acceleration (degree/s^2)');