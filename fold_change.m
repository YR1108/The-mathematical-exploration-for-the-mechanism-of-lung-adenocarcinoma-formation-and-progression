x1=[3.08430000000000,1.09900000000000,1.09970000000000,1.09530000000000,8.06180000000000,1.14150000000000,1.09970000000000,1.00720000000000,5.12530000000000,3.12350000000000,3.03340000000000,0.164000000000000,1.13720000000000,0.1];
x2=[3.16000000000000,1.09900000000000,1.09960000000000,1.09570000000000,7.27560000000000,1.14170000000000,1.09960000000000,1.44000000000000,5.12580000000000,3.12390000000000,3.03350000000000,0.163500000000000,1.13680000000000,0.432900000000000];
x3=[3.18640000000000,1.09900000000000,1.09950000000000,1.09580000000000,6.30620000000000,1.14210000000000,1.09950000000000,1.97320000000000,5.12590000000000,3.12410000000000,3.03360000000000,0.163500000000000,1.13660000000000,0.966100000000000];
p1=x1./x1
p2=x2./x1
p3=x3./x1

figure;
categories = {'EGFR', 'KRAS', 'ALK', 'RET', 'P53', 'P21', 'MET', 'CHEK1', 'AKT', 'CDK2', 'BRAF', 'E2F1', 'ATM', 'ATR'};
bar(categories, p1, 'FaceColor', 'b');
hold on;
bar(categories, p2, 'FaceColor', 'r', 'Offset', [0 0.5 0 0.5 0 0.5 0 0.5 0 0.5 0 0.5 0 0.5]);
hold on;
bar(categories, p3, 'FaceColor', 'g', 'Offset', [0 0.5 0 0.5 0 0.5 0 0.5 0 0.5 0 0.5 0 0.5]);

% 添加图例和标题
legend('Normal', 'premalignant', 'Cancer');
% 显示坐标轴标签
xlabel('genes');
ylabel('fold change');
xticks(1:numel(categories)); % 设置x轴刻度位置
xticklabels(categories); % 设置x轴刻度标签为类别名称
legend('Location', 'NorthWest');

figure('Position', [100, 100, 800, 250]);
x=[1,2,3];
y=[0.001,	0.4329,	0.9661];
% 绘制折线图
plot(x, y,'g','LineWidth',2);
hold on;
scatter(x,y,[],'filled','bo');

% 设置x轴的刻度位置（通常与x变量的值相对应）
xtick_positions = x;

% 设置x轴的刻度标签（你想要显示的变量名称）
% 假设每个x值对应一个名称，例如 'A', 'B', 'C', 'D', 'E'
xtick_labels = {'Normal', 'Premalignant', 'Cancer'};
% 应用x轴的刻度位置和标签
xticks(xtick_positions);
xticklabels(xtick_labels);
title('ATR');

figure('Position', [100, 100, 800, 500]);
b=[1.5, 1.6, 1.7, 1.8, 1.9, 2.0];
n=[1.0095, 0.9653, 0.9471, 0.9352, 0.9252, 0.8854];
p=[0.1383, 0.1784, 0.1888, 0.2017, 0.2153, 0.2251];
c=[0.5796, 0.6132, 0.6291, 0.6352, 0.6358, 0.6822];
% 绘制折线图
% 绘制折线图
plot(b, n,'k-o','LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor','#D95319','DisplayName','Unp');
hold on;
plot(b, p,'k-s','LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor',	'#EDB120','DisplayName','Upn');
hold on;
plot(b, p-0.05,'k-v','LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor', '#77AC30','DisplayName','Upc');
hold on;
plot(b, c,'k-^','LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor','#0072BD','DisplayName','Ucp');
xlabel('repression strength');
ylabel('barrier height');
legend;

figure('Position', [100, 100, 800, 500]);
e=[0.02, 0.03, 0.04, 0.05, 0.06];
n1=[1.6622, 1.3474, 1.1428, 0.9964, 0.8854];
p1=[0.5608, 0.4064, 0.3192, 0.2634, 0.2251];
c1=[1.3699, 1.084, 0.9032, 0.7765, 0.6822];
% 绘制折线图
plot(e, n1,'k-o','LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor','#D95319','DisplayName','U_{N-AAH}');
hold on;
plot(e, p1,'k-s','LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor',	'#EDB120','DisplayName','U_{AAH-N}');
hold on;
plot(e, p1-0.05,'k-v','LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor', '#77AC30','DisplayName','U_{AAH-AIS}');
hold on;
plot(e, c1,'k-^','LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor','#0072BD','DisplayName','U_{AIS-AAH}');
xlabel('noise amplitude');
ylabel('barrier height');
legend;

figure('Position', [100, 100, 800, 250]);
categories = {'EGFR', 'KRAS', 'ALK', 'RET', 'P53', 'P21', 'MET', 'CHEK1', 'AKT', 'CDK2', 'BRAF', 'E2F1', 'ATM', 'ATR'};
Unp=[0.0732, 0.0722, -0.0065, -0.0001, 0.1224, 0.0677, -0.0103, 0.0245, -0.131, 0.0249, -0.0625, 0.0665, -0.0112, -0.3209];
Upn=[-0.0669, -0.0533, 0.0058, -0.0032, -0.1537, -0.0519, -0.0247, 0.0722, 0.0218, 0.0067, 0.0389, -0.016, -0.0183, 0.2168];
bar( Unp, 'FaceColor', '#D95319');
hold on;
bar( Upn, 'FaceColor', '#EDB120');
% 添加图例和标题
legend('\Delta U_{N-AAH}', '\Delta U_{AAH-N}');
% 显示坐标轴标签
xlabel('genes');
ylabel('barrier height change');
xticks(1:numel(categories)); % 设置x轴刻度位置
xticklabels(categories); % 设置x轴刻度标签为类别名称
legend('Location', 'NorthWest');


figure('Position', [100, 100, 800, 250]);
categories = {'EGFR', 'KRAS', 'ALK', 'RET', 'P53', 'P21', 'MET', 'CHEK1', 'AKT', 'CDK2', 'BRAF', 'E2F1', 'ATM', 'ATR'};
Upc=[-0.0669, -0.0533, 0.0058, -0.0032, -0.1537, -0.0519, -0.0247, 0.0722, 0.0218, 0.0067, 0.0389, -0.016, -0.0183, 0.2168];
Ucp=[-0.0939, -0.0822, 0.0082, -0.0035, -0.1216, -0.0756, 0.0271, -0.0696, 0.1519, -0.0371, 0.0664, -0.091, 0.0256, 0.3623];
bar( Upc, 'FaceColor', '#77AC30');
hold on;
bar( Ucp, 'FaceColor', '#0072BD');
% 添加图例和标题
legend('\Delta U_{AAH-AIS}', '\Delta U_{AIS-AAH}');
% 显示坐标轴标签
xlabel('genes');
ylabel('barrier height change');
xticks(1:numel(categories)); % 设置x轴刻度位置
xticklabels(categories); % 设置x轴刻度标签为类别名称
legend('Location', 'NorthWest');


figure('Position', [100, 100, 800, 400]);
categories = {'EGFR→P53', 'KRAS→P53', 'ALK→P53', 'ATM→P53', 'AKT—|P53', 'P21—|P53', 'CDK2—|P53', 'ATR—|P53','ATR→ATR','KRAS→CHEK1', 'ATR→CHEK1', 'AKT—|CHEK1'};
Unp=[0.0913, 0.0849, 0.0894,0.08606,0.0877, 0.1203, 0.1348, 0.1759, -0.1359,0.1165,0.0716,0.1088];
Upn=[-0.0936, -0.0672, -0.061, -0.0505, -0.0867, -0.1101, -0.1207, -0.2696, -0.0838,-0.0872,-0.0443,-0.1111];
colors=['r', 'r', 'r', 'r', 'b', 'b', 'b', 'b', 'r', 'r', 'r', 'b'];
bar( Unp, 'FaceColor', '#D95319');
hold on;
bar( Upn, 'FaceColor', '#EDB120');

% 添加图例和标题
legend('\Delta U_{N-AAH}', '\Delta U_{AAH-N}');
% 显示坐标轴标签
xlabel('regulations');
ylabel('barrier height change');
legend('Location', 'NorthWest');
xticks(1:numel(categories)); % 设置x轴刻度位置
xticklabels(categories); % 设置x轴刻度标签为类别名称

set(gca, 'xticklabelrotation', -90);
% h=xticklabels;
% for i = 1:numel(categories)
%     % 设置每个刻度标签的文本颜色
%     set(h{i,1}, 'Color', colors(i));
% end



figure('Position', [100, 100, 800, 400]);
categories = {'EGFR→P53', 'KRAS→P53', 'ALK→P53', 'ATM→P53', 'AKT—|P53', 'P21—|P53', 'CDK2—|P53', 'ATR—|P53','ATR→ATR','KRAS→CHEK1', 'ATR→CHEK1', 'AKT—|CHEK1'};
Ucp=[-0.0951, -0.0967, -0.1073, -0.0977, -0.0922, -0.1381, -0.1596, -0.1712, 0.1963, -0.1425, -0.086, -0.1173];
Upc=[-0.0936, -0.0672, -0.061, -0.0505, -0.0867, -0.1101, -0.1207, -0.2696, -0.0838,-0.0872, -0.0443, -0.0111];


bar( Upc, 'FaceColor', '#77AC30');
hold on;
bar( Ucp, 'FaceColor', '#0072BD');

% 添加图例和标题
legend('\Delta U_{AAH-AIS}', '\Delta U_{AIS-AAH}');
% 显示坐标轴标签
xlabel('regulations');
ylabel('barrier height change');
xticks(1:numel(categories)); % 设置x轴刻度位置
xticklabels(categories); % 设置x轴刻度标签为类别名称
legend('Location', 'NorthWest');
set(gca, 'xticklabelrotation', -90);

%% basin 热图
figure('Position', [100, 100, 1200, 400]);
% 绘制热度图
% imagesc(newnetwork1); % 或者使用 pimagesc(data) 来获得更好的性能
% 添加颜色条
% colorbar;
Z = zeros(11,31);
imagesc(Z);
% 添加轴标签和标题
xlabel('Activation  constant');
ylabel('Repression constant');

% 设置x轴和y轴的刻度位置和标签
x = 1:31; % 设置x轴的刻度位置
y = 1:11; % 设置y轴的刻度位置
x_ticks = 0.5:0.05:2; 
y_ticks = 1:0.1:2;
% 为x轴和y轴添加刻度
xticks(x);
yticks(y);
% 为x轴和y轴添加刻度标签
xticklabels(arrayfun(@num2str, x_ticks, 'UniformOutput', false)); % 将数字转换为字符串作为标签
yticklabels(arrayfun(@num2str, y_ticks, 'UniformOutput', false));

% 可选：设置颜色映射
% colormap('summer'); % 使用'hot'颜色映射，也可以选择其他颜色映射

%% basin 热图
figure('Position', [50, 20, 400, 1000]);
% 绘制热度图
imagesc(newnetwork1'); % 或者使用 pimagesc(data) 来获得更好的性能

% 添加颜色条
colorbar;

% 添加轴标签和标题
ylabel('Activation  constant');
xlabel('Repression constant');

% 设置x轴和y轴的刻度位置和标签
y = 1:31; % 设置x轴的刻度位置
x = 1:11; % 设置y轴的刻度位置

% 为x轴和y轴添加刻度
yticks(y);
xticks(x);

% 为x轴和y轴添加刻度标签
xticklabels(arrayfun(@num2str, y_ticks, 'UniformOutput', false)); % 将数字转换为字符串作为标签
yticklabels(arrayfun(@num2str, x_ticks, 'UniformOutput', false));

% 可选：设置颜色映射
colormap('summer'); % 使用'hot'颜色映射，也可以选择其他颜色映射

