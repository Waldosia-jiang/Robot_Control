#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

// 参数定义
const double k = 0.1;       // 前视增益，影响前视距离随速度的动态调整
const double Lfc = 2.0;     // 基础前视距离 [m]，确保低速时有足够的前视点
const double Kp = 1.0;      // 速度比例增益，用于比例控制加速度
const double dt = 0.1;      // 时间步长 [s]，模拟的时间分辨率
const double WB = 2.9;      // 车辆轴距 [m]，车辆前后轮之间的距离
const double max_lat_acc = 2.0; // 最大横向加速度 [m/s^2]，限制转弯时的速度
const double initial_target_speed = 10.0 / 3.6; // 初始目标速度 [m/s]，转换为 10 km/h
const bool show_animation = true; // 是否显示动画，用于实时可视化

// 车辆状态类，管理车辆的位置、速度和航向
class State {
public:
    double x, y, yaw, v, rear_x, rear_y; // 位置(x, y)，偏航角(yaw)，速度(v)，后轮位置(rear_x, rear_y)

    // 构造函数，初始化车辆状态
    State(double x_ = 0.0, double y_ = 0.0, double yaw_ = 0.0, double v_ = 0.0)
        : x(x_), y(y_), yaw(yaw_), v(v_) {
        updateRear(); // 初始化后轮位置
    }

    // 更新车辆状态，根据加速度和转向角
    void update(double a, double delta) {
        x += v * std::cos(yaw) * dt;       // 更新x坐标，基于速度和航向角
        y += v * std::sin(yaw) * dt;       // 更新y坐标，基于速度和航向角
        yaw += v / WB * std::tan(delta) * dt; // 更新偏航角，基于单车模型的角速度
        v += a * dt;                       // 更新速度，基于加速度
        updateRear();                      // 更新后轮位置
    }

    // 计算后轮与某点的欧几里得距离
    double calcDistance(double point_x, double point_y) const {
        double dx = rear_x - point_x; // x方向差值
        double dy = rear_y - point_y; // y方向差值
        return std::hypot(dx, dy);    // 返回两点间的距离
    }

private:
    // 更新后轮位置，基于车辆中心位置和航向角
    void updateRear() {
        rear_x = x - (WB / 2) * std::cos(yaw); // 后轮x坐标
        rear_y = y - (WB / 2) * std::sin(yaw); // 后轮y坐标
    }
};

// 存储车辆状态历史，用于记录轨迹和分析
class States {
public:
    std::vector<double> x, y, yaw, v, t; // 存储x, y，偏航角，速度，时间的历史记录

    // 添加当前状态到历史记录
    void append(double time, const State& state) {
        x.push_back(state.x);
        y.push_back(state.y);
        yaw.push_back(state.yaw);
        v.push_back(state.v);
        t.push_back(time);
    }
};

// PID速度控制，仅使用比例控制（P）计算加速度
double proportionalControl(double target, double current) {
    return Kp * (target - current); // 加速度 = 比例增益 * (目标速度 - 当前速度)
}

// 计算基于转弯半径的最大速度，限制转弯时的速度
double calcMaxSpeed(double delta) {
    if (std::fabs(delta) < 1e-5) return initial_target_speed; // 如果转向角很小（接近直线），返回初始目标速度
    double R = WB / std::tan(delta); // 计算转弯半径，基于单车模型
    return std::sqrt(max_lat_acc * std::fabs(R)); // 返回最大速度，v_max = sqrt(a_lat_max * |R|)
}

// 目标路径类，管理路径点和目标点搜索
class TargetCourse {
public:
    std::vector<double> cx, cy; // 目标路径的x, y坐标
    int old_nearest_point_index; // 上一次最近点的索引，用于优化搜索

    // 构造函数，初始化路径坐标
    TargetCourse(const std::vector<double>& cx_, const std::vector<double>& cy_)
        : cx(cx_), cy(cy_), old_nearest_point_index(-1) {}

    // 搜索目标点索引和前视距离
    std::pair<int, double> searchTargetIndex(const State& state) {
        int ind;
        if (old_nearest_point_index == -1) { // 第一次搜索最近点
            std::vector<double> d; // 存储所有点与后轮的距离
            for (size_t i = 0; i < cx.size(); ++i) {
                double dx = state.rear_x - cx[i];
                double dy = state.rear_y - cy[i];
                d.push_back(std::hypot(dx, dy));
            }
            ind = std::min_element(d.begin(), d.end()) - d.begin(); // 找到最近点的索引
            old_nearest_point_index = ind; // 记录最近点
        } else { // 从上一次最近点开始搜索
            ind = old_nearest_point_index;
            double distance_this_index = state.calcDistance(cx[ind], cy[ind]); // 当前点的距离
            while (true) {
                if (ind + 1 >= static_cast<int>(cx.size())) break; // 超出路径范围则停止
                double distance_next_index = state.calcDistance(cx[ind + 1], cy[ind + 1]); // 下一点的距离
                if (distance_this_index < distance_next_index) break; // 如果当前点更近，则停止
                ind++; // 移动到下一点
                distance_this_index = distance_next_index;
            }
            old_nearest_point_index = ind; // 更新最近点索引
        }

        double Lf = k * state.v + Lfc; // 计算动态前视距离，随速度增加而变长
        while (Lf > state.calcDistance(cx[ind], cy[ind])) { // 寻找前视点
            if (ind + 1 >= static_cast<int>(cx.size())) break; // 超出路径范围则停止
            ind++; // 移动到下一点
        }
        return {ind, Lf}; // 返回目标点索引和前视距离
    }
};

// Pure Pursuit转向控制，计算转向角和目标点索引
std::pair<double, int> purePursuitSteerControl(const State& state, TargetCourse& trajectory, int pind) {
    std::pair<int, double> target_info = trajectory.searchTargetIndex(state); // 获取目标点信息
    int ind = target_info.first;    // 目标点索引
    double Lf = target_info.second; // 前视距离

    if (pind >= ind) ind = pind; // 确保目标点索引不会回退

    // 获取目标点坐标，若超出范围则使用路径终点
    double tx = (ind < static_cast<int>(trajectory.cx.size())) ? trajectory.cx[ind] : trajectory.cx.back();
    double ty = (ind < static_cast<int>(trajectory.cy.size())) ? trajectory.cy[ind] : trajectory.cy.back();

    double alpha = std::atan2(ty - state.rear_y, tx - state.rear_x) - state.yaw; // 计算目标点与当前航向的夹角
    double delta = std::atan2(2.0 * WB * std::sin(alpha) / Lf, 1.0); // 计算转向角，基于纯追踪几何公式

    return {delta, ind}; // 返回转向角和目标点索引
}

// 绘制箭头，表示车辆位置和方向
void plotArrow(double x, double y, double yaw, double length = 1.0, double width = 0.5) {
    std::vector<double> x_vec = {x, x + length * std::cos(yaw)}; // 箭头起止点的x坐标
    std::vector<double> y_vec = {y, y + length * std::sin(yaw)}; // 箭头起止点的y坐标
    plt::plot(x_vec, y_vec, "r-"); // 绘制箭身（红线）
    plt::plot({x}, {y}, "ro");     // 绘制起点（红点）
}

// 主模拟函数，执行路径跟踪模拟
void mainSimulation() {
    std::vector<double> cx, cy; // 目标路径坐标
    for (double i = 0; i < 50; i += 0.5) { // 生成正弦路径
        cx.push_back(i);
        cy.push_back(std::sin(i / 5.0) * i / 2.0);
    }

    State state(-0.0, -3.0, 0.0, 0.0); // 初始化车辆状态：起点在(0, -3)，航向0，速度0
    States states; // 创建状态记录对象
    states.append(0.0, state); // 记录初始状态

    TargetCourse target_course(cx, cy); // 创建目标路径对象
    int target_ind = target_course.searchTargetIndex(state).first; // 初始目标点索引
    int lastIndex = cx.size() - 1;    // 路径终点索引
    double time = 0.0;                // 当前时间
    const double T = 100.0;           // 最大模拟时间

    while (T >= time && lastIndex > target_ind) { // 主循环：未超时且未到达终点
        // 计算转向角和目标点
        std::pair<double, int> steer_info = purePursuitSteerControl(state, target_course, target_ind);
        double di = steer_info.first;   // 转向角
        target_ind = steer_info.second; // 更新目标点索引

        // 根据转弯半径计算最大速度，并限制目标速度
        double max_speed = calcMaxSpeed(di);
        double target_speed = std::min(initial_target_speed, max_speed);

        // 计算加速度
        double ai = proportionalControl(target_speed, state.v);

        // 更新车辆状态
        state.update(ai, di);
        time += dt;           // 时间递增
        states.append(time, state); // 记录当前状态

        if (show_animation) { // 实时动画显示
            plt::cla();       // 清空当前画布
            plotArrow(state.x, state.y, state.yaw); // 绘制车辆方向箭头
            std::map<std::string, std::string> keywords_course = {{"label", "course"}, {"color", "red"}, {"linestyle", "-"}};
            plt::plot(cx, cy, keywords_course); // 绘制目标路径（红线）
            std::map<std::string, std::string> keywords_traj = {{"label", "trajectory"}, {"color", "blue"}, {"linestyle", "-"}};
            plt::plot(states.x, states.y, keywords_traj); // 绘制实际轨迹（蓝线）
            std::map<std::string, std::string> keywords_target = {{"label", "target"}, {"marker", "x"}, {"color", "green"}};
            plt::plot({cx[target_ind]}, {cy[target_ind]}, keywords_target); // 绘制目标点（绿色x）
            plt::axis("equal"); // 设置坐标轴等比例
            plt::grid(true);    // 显示网格
            plt::title("速度[km/h]: " + std::to_string(state.v * 3.6).substr(0, 4)); // 显示当前速度
            plt::pause(0.001);  // 暂停以显示动画帧
        }
    }

    if (show_animation) { // 显示最终结果
        plt::cla();       // 清空画布
        std::map<std::string, std::string> keywords_course = {{"label", "course"}, {"color", "red"}, {"marker", "."}};
        plt::plot(cx, cy, keywords_course); // 绘制目标路径（红点）
        std::map<std::string, std::string> keywords_traj = {{"label", "trajectory"}, {"color", "blue"}, {"linestyle", "-"}};
        plt::plot(states.x, states.y, keywords_traj); // 绘制实际轨迹（蓝线）
        plt::legend();    // 显示图例
        plt::xlabel("x[m]"); // x轴标签
        plt::ylabel("y[m]"); // y轴标签
        plt::axis("equal");  // 设置坐标轴等比例
        plt::grid(true);     // 显示网格

        plt::figure();       // 创建新画布用于速度曲线
        std::map<std::string, std::string> keywords_speed = {{"color", "red"}, {"linestyle", "-"}};
        std::vector<double> speed_kmh(states.v.size()); // 速度转换为km/h
        for (size_t i = 0; i < states.v.size(); ++i) speed_kmh[i] = states.v[i] * 3.6;
        plt::plot(states.t, speed_kmh, keywords_speed); // 绘制速度-时间曲线
        plt::xlabel("时间[s]"); // x轴标签
        plt::ylabel("速度[km/h]"); // y轴标签
        plt::grid(true);     // 显示网格
        plt::show();         // 显示所有图形
    }
}

// 主函数，程序入口
int main() {
    std::cout << "Pure pursuit路径跟踪模拟开始" << std::endl; // 打印开始信息
    mainSimulation(); // 运行主模拟函数
    return 0;         // 程序正常退出
}
