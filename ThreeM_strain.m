clc; clear; close all;
bear_data = readtable("轴系参数表_pilot.xlsx",'VariableNamingRule','preserve');

number_sections = size(bear_data,1);%设置截面数量

%设置矩阵大小
A11 = zeros(number_sections,number_sections);A12 = zeros(number_sections,number_sections);
A21 = zeros(2,number_sections);A22 = zeros(2,number_sections);
A31 = zeros(number_sections,number_sections);A32 = zeros(number_sections,number_sections);
C = zeros(1,2*number_sections);
%设置L，E，I，Q；
L = zeros(1,number_sections);
E = zeros(1,number_sections);
I = zeros(1,number_sections);
Q = zeros(1,number_sections);
Z = zeros(1,number_sections);
Position_Bearing = zeros(1,number_sections);
PointLoad = zeros(1,number_sections);
for i = 1:number_sections
    L(i) = bear_data.("单元长度mm")(i);%单位mm
    E(i) = bear_data.("弹性模量GPa")(i);%弹性模量单位KN/mm2=Gpa
    I(i) = bear_data.("右面外径mm")(i)^4*pi/64;%截面惯性矩mm4
    Q(i) = 0.25*bear_data.("材料密度Kg/m^3")(i)*bear_data.("右面外径mm")(i)^2*pi*9.8*10^(-12);%均布载荷KN/mm
    PointLoad(i) = bear_data.("集中载荷kN")(i);
    PointLoad(i) = bear_data.("集中载荷kN")(i);
    Displacement_Bearing(i) = bear_data.("轴承变位")(i);
    Position_Bearing(i) = bear_data.("轴承位置（1/0）")(i);
end

%一段自由端
%填充矩阵A11,A12
% 合并A21，22到A11，12
for i = 1:number_sections
    if i == 1
        A11(i,i) = 1;
    elseif i == number_sections
        A11(i,i) = 1;
    else
        A11(i,i-1) = L(i-1)/E(i-1)/I(i-1);
        A11(i,i+1) = L(i)/E(i)/I(i);
        A11(i,i) = 2*(A11(i,i-1)+A11(i,i+1));

        A12(i,i-1) = -6/L(i-1);
        A12(i,i+1) = -6/L(i);
        A12(i,i) = -A12(i,i-1)-A12(i,i+1);
    end
end

%填充矩阵A31，A32
for i = 1:number_sections
    if i == 1
        A31(i,i) = -1/L(i);
        A31(i,i+1) = 1/L(i);

    elseif i == number_sections
        A31(number_sections,number_sections-1) =1/L(i-1);
        A31(number_sections,number_sections) =-1/L(i-1);
        %elseif i==number_sections-1

    elseif i>1

        A31(i,i-1) = 1/L(i-1);
        A31(i,i+1) = 1/L(i);
        A31(i,i) = -A31(i,i-1)-A31(i,i+1);

        if Position_Bearing(i) == 1
            A31(i,i-1) = 0;
            A31(i,i+1) = 0;
            A31(i,i) = 0;
        end
    end
    A32(i,i) = Position_Bearing(i);
end

%% 此C矩阵第一个为合并A11，A21，A12，A22后的C
for i = 1:number_sections
    if i == 1
        C(i) = 0;
    elseif i == number_sections
        C(i) = 0;
    else
        C(i) = -0.25*(Q(i)*L(i)^3/(E(i)*I(i))+Q(i-1)*L(i-1)^3/(E(i-1)*I(i-1)));
    end
end
% C(number_sections) = -0.25*(Q(i-1)*L(i-1)^3/(E(i11)*I(i-1)));%若为固定端则关闭注释

for i = number_sections+1:2*number_sections
    j=i-number_sections;
    if Position_Bearing(j)==0
        if j == 1
            C(i) = -0.5*L(j)*Q(j)-PointLoad(j);
        elseif j== number_sections
            C(i) = -0.5*L(j-1)*Q(j-1)-PointLoad(j);
        else
            C(i) = -0.5*(L(j-1)*Q(j-1)+L(j)*Q(j))-PointLoad(j);
        end
    else
        C(i) = Displacement_Bearing(j);
    end
end
%计算矩阵A*X=C
A = [A11 A12;A31 A32];
X = inv(A)*C';%弯矩KN·mm=n*m  %挠度：弯曲变形时横截面形心沿与轴线垂直方向的线位移
F = nan(number_sections,1);
for i = 2:number_sections-1
    if bear_data.("轴承位置（1/0）")(i)==1
        F(i) = (X(i-1)-X(i))/L(i-1)+(X(i+1)-X(i))/L(i)+Q(i-1).*L(i-1)/2+Q(i).*L(i)/2+PointLoad(i);
        continue
    end
end

F_list = F(~isnan(F));
F_sum = sum(F_list);
for i = 1:number_sections
    y1(i) = X(i);
    z1(i) = X(i+38);
    t1(i) = bear_data.("节点坐标mm")(i);
end
figure(1)
plot(t1,F,'*')
legend('力/KN')
% plot(t1,F,'b--o')
% legend('力（KN）')
% xlabel('l');
% ylabel('y');
% title('三弯矩方程')
figure(2)
X_size = size(X,1)/2+1;
plot(t1,X(X_size:end,1),'*')
hold on 
plot(t1,X(X_size:end,1))
hold off

