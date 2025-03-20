#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

// 参数定义（全局常量）
const double k = 1.0;         // Stanley 控制增益，调节横向误差对转向的影响
const double Kp = 1.0;        // PID 速度控制比例增益
const double dt = 0.1;        // 时间步长 [s]，模拟的采样周期
const double L = 2.9;         // 车辆轴距 [m]，用于运动学模型
const double max_steer = M_PI / 6.0; // 最大转向角 [rad]，约 30 度
const double approach_threshold = 0.5; // 靠近目标点的距离阈值 [m]，决定何时切换到跟踪
const double max_lat_acc = 3.0; // 最大横向加速度 [m/s^2]，用于转弯速度限制
const double initial_target_speed = 10.0 / 3.6; // 初始目标速度 [m/s]，10 km/h
const bool show_animation = true; // 是否显示实时动画

#ifndef M_PI
#define M_PI 3.14159265358979323846 // 定义 PI，如果未预定义
#endif

// 角度规范化函数：将角度限制在 [-π, π] 范围内
double normalize_angle(double angle) {
    while (angle > M_PI) angle -= 2.0 * M_PI;
    while (angle < -M_PI) angle += 2.0 * M_PI;
    return angle;
}

// 车辆状态类：封装位置、航向和速度
class State {
public:
    double x, y;   // 车辆位置坐标 [m]
    double yaw;    // 航向角 [rad]
    double v;      // 速度 [m/s]

    // 构造函数，初始化状态
    State(double x_ = 0.0, double y_ = 0.0, double yaw_ = 0.0, double v_ = 0.0)
        : x(x_), y(y_), yaw(yaw_), v(v_) {}

    // 更新车辆状态基于运动学模型
    void update(double acceleration, double delta) {
        delta = std::max(-max_steer, std::min(max_steer, delta)); // 限制转向角
        x += v * std::cos(yaw) * dt;  // 更新 x 坐标
        y += v * std::sin(yaw) * dt;  // 更新 y 坐标
        yaw += v / L * std::tan(delta) * dt; // 更新航向角（自行车模型）
        yaw = normalize_angle(yaw);   // 规范化航向角
        v += acceleration * dt;       // 更新速度
    };

    // 计算到指定点的欧几里得距离
    double distance_to(double px, double py) const {
        double dx = x - px;
        double dy = y - py;
        return std::hypot(dx, dy);
    }
};

// PID 速度控制：简单比例控制，计算加速度
double pid_control(double target, double current) {
    return Kp * (target - current); // 比例控制输出加速度
}

// 计算基于转弯半径的最大速度，限制转弯时的速度
double calc_max_speed(double delta) {
    if (std::fabs(delta) < 1e-5) return initial_target_speed; // 如果转向角很小（接近直线），返回初始目标速度
    double R = L / std::tan(delta); // 计算转弯半径，基于单车模型
    return std::sqrt(max_lat_acc * std::fabs(R)); // 返回最大速度，v_max = sqrt(a_lat_max * |R|)
}

// 计算目标点索引和前轴横向误差
std::pair<int, double> calc_target_index(const State& state, const std::vector<double>& cx, const std::vector<double>& cy) {
    // 计算前轴位置（基于车辆轴距和航向）
    double fx = state.x + L * std::cos(state.yaw);
    double fy = state.y + L * std::sin(state.yaw);

    // 计算前轴到路径所有点的距离
    std::vector<double> d;
    for (size_t i = 0; i < cx.size(); ++i) {
        double dx = fx - cx[i];
        double dy = fy - cy[i];
        d.push_back(std::hypot(dx, dy));
    }
    // 找到最近点的索引
    int target_idx = std::distance(d.begin(), std::min_element(d.begin(), d.end()));

    // 计算横向误差（前轴到最近点的垂直距离）
    std::vector<double> front_axle_vec = {-std::cos(state.yaw + M_PI / 2), -std::sin(state.yaw + M_PI / 2)}; // 前轴法向量
    double dx = fx - cx[target_idx];
    double dy = fy - cy[target_idx];
    double error_front_axle = dx * front_axle_vec[0] + dy * front_axle_vec[1]; // 点积计算横向误差

    return {target_idx, error_front_axle};
}

// Stanley 转向控制：计算转向角和新的目标点索引
std::pair<double, int> stanley_control(const State& state, const std::vector<double>& cx, 
                                       const std::vector<double>& cy, const std::vector<double>& cyaw, int last_target_idx) {
    // 获取当前目标点索引和横向误差
    auto [current_target_idx, error_front_axle] = calc_target_index(state, cx, cy);

    // 确保目标点不回退（单向前进）
    if (last_target_idx >= current_target_idx) {
        current_target_idx = last_target_idx;
    }

    // 计算航向误差（路径方向与车辆方向的差）
    double theta_e = normalize_angle(cyaw[current_target_idx] - state.yaw);
    // 计算横向误差校正项（基于 Stanley 公式）
    double theta_d = std::atan2(k * error_front_axle, state.v + 1e-6); // 加小量避免除零
    // 总转向角
    double delta = theta_e + theta_d;

    return {delta, current_target_idx};
}

// 简单转向控制：用于初始靠近阶段，基于几何转向
double simple_steering_control(const State& state, double target_x, double target_y) {
    double target_yaw = std::atan2(target_y - state.y, target_x - state.x); // 目标方向
    double yaw_error = normalize_angle(target_yaw - state.yaw); // 航向误差
    return std::atan2(2.0 * L * std::sin(yaw_error), 1.0); // 简单几何转向公式
}

// 生成参考路径：与 Pure Pursuit 相同的正弦路径
void generate_reference_path(std::vector<double>& cx, std::vector<double>& cy, std::vector<double>& cyaw) {
    // 生成与 Pure Pursuit 相同的正弦路径
    for (double x = 0; x < 50; x += 0.5) {
        cx.push_back(x);
        cy.push_back(std::sin(x / 5.0) * x / 2.0); // y = sin(x/5) * x/2
    }

    // 计算路径方向（航向角）
    cyaw.resize(cx.size());
    for (size_t i = 0; i < cx.size() - 1; ++i) {
        double dx = cx[i + 1] - cx[i];
        double dy = cy[i + 1] - cy[i];
        cyaw[i] = std::atan2(dy, dx); // 两点间的切线方向
    }
    cyaw.back() = cyaw[cyaw.size() - 2]; // 最后一个点方向与前一个相同

    // 平滑航向角（简单平均）
    std::vector<double> cyaw_smooth(cx.size());
    cyaw_smooth[0] = cyaw[0];
    cyaw_smooth.back() = cyaw.back();
    for (size_t i = 1; i < cx.size() - 1; ++i) {
        cyaw_smooth[i] = (cyaw[i - 1] + cyaw[i] + cyaw[i + 1]) / 3.0; // 三点平均平滑
    }
    cyaw = cyaw_smooth;
}

// 绘制箭头：显示车辆位置和航向
void plot_arrow(double x, double y, double yaw, double length = 1.0) {
    std::vector<double> x_vec = {x, x + length * std::cos(yaw)}; // 箭头起止点的 x 坐标
    std::vector<double> y_vec = {y, y + length * std::sin(yaw)}; // 箭头起止点的 y 坐标
    plt::plot(x_vec, y_vec, "r-"); // 绘制箭身（红线）
    plt::plot({x}, {y}, "ro");     // 绘制起点（红点）
}

// 主模拟函数：执行路径跟踪模拟
void main_simulation() {
    // 生成参考路径
    std::vector<double> cx, cy, cyaw;
    generate_reference_path(cx, cy, cyaw);

    State state(-0.0, -3.0, 0.0, 0.0); // 初始状态：位置 (-0, -3)，航向 0，速度 0

    int target_idx = calc_target_index(state, cx, cy).first; // 初始目标点索引
    int last_idx = cx.size() - 1; // 路径终点索引
    bool is_tracking = false; // 是否进入跟踪模式
    double time = 0.0; // 当前模拟时间
    std::vector<double> x = {state.x}; // 记录 x 坐标轨迹
    std::vector<double> y = {state.y}; // 记录 y 坐标轨迹
    std::vector<double> yaw = {state.yaw}; // 记录航向轨迹
    std::vector<double> v = {state.v}; // 记录速度轨迹
    std::vector<double> t = {0.0}; // 记录时间轨迹

    std::cout << "Initial target index: " << target_idx << ", coordinates: (" << cx[target_idx] << ", " << cy[target_idx] << ")" << std::endl;

    // 主模拟循环
    while (time < 100.0 && (!is_tracking || last_idx > target_idx)) {
        // 计算转向角
        double delta;
        if (!is_tracking) {
            // 靠近阶段：动态更新目标点并移动到最近点
            target_idx = calc_target_index(state, cx, cy).first;
            delta = simple_steering_control(state, cx[target_idx], cy[target_idx]); // 计算简单转向角
        } else {
            // Stanley 跟踪阶段：使用 Stanley 控制跟踪路径
            auto [di, new_target_idx] = stanley_control(state, cx, cy, cyaw, target_idx);
            delta = di;
            target_idx = new_target_idx; // 更新目标点索引
        }

        // 根据转弯半径计算最大速度，并限制目标速度
        double max_speed = calc_max_speed(delta);
        double target_speed = std::min(initial_target_speed, max_speed);

        // 计算加速度
        double ai = pid_control(target_speed, state.v);

        // 更新车辆状态
        state.update(ai, delta);

        // 检查是否进入跟踪模式
        if (!is_tracking) {
            double dist_to_target = state.distance_to(cx[target_idx], cy[target_idx]); // 计算到目标点的距离
            if (dist_to_target < approach_threshold) {
                is_tracking = true; // 距离足够近，切换到跟踪模式
                std::cout << "Reached target " << target_idx << " (" << cx[target_idx] << ", " << cy[target_idx] 
                          << "), starting path tracking" << std::endl;
            }
        }

        time += dt; // 时间递增
        x.push_back(state.x); // 记录轨迹
        y.push_back(state.y);
        yaw.push_back(state.yaw);
        v.push_back(state.v);
        t.push_back(time);

        // 实时动画显示
        if (show_animation) {
            plt::cla(); // 清空当前画布
            std::map<std::string, std::string> keywords_course = {{"label", "course"}, {"color", "red"}, {"marker", "."}};
            plt::plot(cx, cy, keywords_course); // 绘制参考路径（红点）
            std::map<std::string, std::string> keywords_traj = {{"label", "trajectory"}, {"color", "blue"}, {"linestyle", "-"}};
            plt::plot(x, y, keywords_traj); // 绘制车辆轨迹（蓝线）
            std::map<std::string, std::string> keywords_target = {{"label", "target"}, {"marker", "x"}, {"color", "green"}};
            plt::plot({cx[target_idx]}, {cy[target_idx]}, keywords_target); // 绘制目标点（绿叉）
            plot_arrow(state.x, state.y, state.yaw); // 绘制航向箭头
            plt::axis("equal"); // 设置坐标轴等比例
            plt::grid(true); // 显示网格
            plt::xlabel("X [m]"); // x 轴标签
            plt::ylabel("Y [m]"); // y 轴标签
            plt::title("Speed [km/h]: " + std::to_string(state.v * 3.6).substr(0, 4) + 
                       (is_tracking ? " (Tracking)" : " (Approaching)")); // 标题显示速度和状态
            plt::pause(0.001); // 暂停以刷新动画
        }
    }

    // 检查是否到达终点
    if (last_idx <= target_idx) {
        std::cout << "Reached the end of the path!" << std::endl;
    } else {
        std::cerr << "Failed to reach the end of the path" << std::endl;
    }

    // 最终结果可视化
    if (show_animation) {
        plt::cla();
        std::map<std::string, std::string> keywords_course = {{"label", "course"}, {"color", "red"}, {"marker", "."}};
        plt::plot(cx, cy, keywords_course); // 绘制参考路径
        std::map<std::string, std::string> keywords_traj = {{"label", "trajectory"}, {"color", "blue"}, {"linestyle", "-"}};
        plt::plot(x, y, keywords_traj); // 绘制完整轨迹
        plt::legend(); // 显示图例
        plt::xlabel("X [m]");
        plt::ylabel("Y [m]");
        plt::axis("equal");
        plt::grid(true);

        plt::figure(); // 新建速度图窗口
        std::map<std::string, std::string> keywords_speed = {{"color", "red"}, {"linestyle", "-"}};
        std::vector<double> speed_kmh(v.size());
        for (size_t i = 0; i < v.size(); ++i) {
            speed_kmh[i] = v[i] * 3.6; // 转换为 km/h
        }
        plt::plot(t, speed_kmh, keywords_speed); // 绘制速度曲线
        plt::xlabel("Time [s]");
        plt::ylabel("Speed [km/h]");
        plt::grid(true);
        plt::show(); // 显示所有图形
    }
}

// 主函数：程序入口
int main() {
    std::cout << "Stanley path tracking simulation started" << std::endl;
    main_simulation(); // 运行模拟
    return 0;
}