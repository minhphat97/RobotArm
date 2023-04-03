data = csvread('joint_space_position_segment.csv');
time = data(:,1);

joint_position_1 = data(:,2);
figure(1);
plot(time, joint_position_1, '-o');
title('Position of Joint 1, 2 & 4');
xlabel('Time (s)');
ylabel('Position (degree)');
hold on;

joint_position_2 = data(:,3);
plot(time, joint_position_2, '-ro');

joint_position_4 = data(:,5);
plot(time, joint_position_4, '-go');
legend("position of joint 1", "position of joint 2", "position of joint 4");
hold off;

joint_position_3 = data(:,4);
figure(2);
plot(time, joint_position_3, '-o');
title('Position of Joint 3');
xlabel('Time (s)');
ylabel('Position (mm)');

x_data = data(:,6);
y_data = data(:,7);
figure(3);
plot(x_data, y_data, '-o');
title('X vs Y');
xlim([-400 400]);
ylim([-400 400]);
xlabel('X (mm)');
ylabel('Y (mm)');