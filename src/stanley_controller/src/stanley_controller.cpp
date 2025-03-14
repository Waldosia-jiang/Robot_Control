#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include <matplotlibcpp.h>

namespace plt = matplotlibcpp;

// 参数定义
const double k = 0.5;         // Stanley控制增益
const double Kp = 1.0;        // 速度比例增益
const double dt = 0.1;        // 时间步长 [s]
const double L = 2.9;         // 车辆轴距 [m]
const double max_steer = M_PI / 6.0; // 最大转向角 [rad]，30度
const double approach_threshold = 0.5; // 靠近目标点的距离阈值 [m]
const bool show_animation = true; // 是否显示动画

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// 角度规范化函数
double normalize_angle(double angle) {
    while (angle > M_PI) angle -= 2.0 * M_PI;
    while (angle < -M_PI) angle += 2.0 * M_PI;
    return angle;
}

// 车辆状态类
class State {
public:
    double x, y, yaw, v;

    State(double x_ = 0.0, double y_ = 0.0, double yaw_ = 0.0, double v_ = 0.0)
        : x(x_), y(y_), yaw(yaw_), v(v_) {}

    void update(double acceleration, double delta) {
        delta = std::max(-max_steer, std::min(max_steer, delta));
        x += v * std::cos(yaw) * dt;
        y += v * std::sin(yaw) * dt;
        yaw += v / L * std::tan(delta) * dt;
        yaw = normalize_angle(yaw);
        v += acceleration * dt;
    }

    double distance_to(double px, double py) const {
        double dx = x - px;
        double dy = y - py;
        return std::hypot(dx, dy);
    }
};

// PID速度控制
double pid_control(double target, double current) {
    return Kp * (target - current);
}

// 计算目标点索引和前轴误差
std::pair<int, double> calc_target_index(const State& state, const std::vector<double>& cx, const std::vector<double>& cy) {
    double fx = state.x + L * std::cos(state.yaw);
    double fy = state.y + L * std::sin(state.yaw);

    std::vector<double> d;
    for (size_t i = 0; i < cx.size(); ++i) {
        double dx = fx - cx[i];
        double dy = fy - cy[i];
        d.push_back(std::hypot(dx, dy));
    }
    int target_idx = std::distance(d.begin(), std::min_element(d.begin(), d.end()));

    std::vector<double> front_axle_vec = {-std::cos(state.yaw + M_PI / 2), -std::sin(state.yaw + M_PI / 2)};
    double dx = fx - cx[target_idx];
    double dy = fy - cy[target_idx];
    double error_front_axle = dx * front_axle_vec[0] + dy * front_axle_vec[1];

    return {target_idx, error_front_axle};
}

// Stanley转向控制
std::pair<double, int> stanley_control(const State& state, const std::vector<double>& cx, 
                                       const std::vector<double>& cy, const std::vector<double>& cyaw, int last_target_idx) {
    auto [current_target_idx, error_front_axle] = calc_target_index(state, cx, cy);

    if (last_target_idx >= current_target_idx) {
        current_target_idx = last_target_idx;
    }

    double theta_e = normalize_angle(cyaw[current_target_idx] - state.yaw);
    double theta_d = std::atan2(k * error_front_axle, state.v + 1e-6);
    double delta = theta_e + theta_d;

    return {delta, current_target_idx};
}

// 简单转向控制（用于靠近目标点）
double simple_steering_control(const State& state, double target_x, double target_y) {
    double target_yaw = std::atan2(target_y - state.y, target_x - state.x);
    double yaw_error = normalize_angle(target_yaw - state.yaw);
    return std::atan2(2.0 * L * std::sin(yaw_error), 1.0);
}

// 生成近似直线的缓和路径
void generate_wiggly_straight_path(std::vector<double>& cx, std::vector<double>& cy, std::vector<double>& cyaw) {
    const double total_length = 50.0;
    const double ds = 0.05;
    const double slope = 0.1;
    const double wiggle_amplitude = 1.0;
    const double wiggle_freq = 0.2;

    for (double t = 0; t <= total_length; t += ds) {
        double x = t;
        double y = slope * t + wiggle_amplitude * std::sin(wiggle_freq * t);
        cx.push_back(x);
        cy.push_back(y);
    }

    cyaw.resize(cx.size());
    for (size_t i = 0; i < cx.size() - 1; ++i) {
        double dx = cx[i + 1] - cx[i];
        double dy = cy[i + 1] - cy[i];
        cyaw[i] = std::atan2(dy, dx);
    }
    cyaw.back() = cyaw[cyaw.size() - 2];

    std::vector<double> cyaw_smooth(cx.size());
    cyaw_smooth[0] = cyaw[0];
    cyaw_smooth.back() = cyaw.back();
    for (size_t i = 1; i < cx.size() - 1; ++i) {
        cyaw_smooth[i] = (cyaw[i - 1] + cyaw[i] + cyaw[i + 1]) / 3.0;
    }
    cyaw = cyaw_smooth;

    std::cout << "First 5 path points:" << std::endl;
    for (size_t i = 0; i < std::min<size_t>(5, cx.size()); ++i) {
        std::cout << "x: " << cx[i] << ", y: " << cy[i] << ", yaw: " << cyaw[i] << std::endl;
    }
}

// 主模拟函数
void main_simulation() {
    // 生成路径
    std::vector<double> cx, cy, cyaw;
    generate_wiggly_straight_path(cx, cy, cyaw);

    double target_speed = 30.0 / 3.6;
    double max_simulation_time = 100.0;

    State state(-0.0, -3.0, 0.0, 0.0);

    int target_idx = calc_target_index(state, cx, cy).first;
    int last_idx = cx.size() - 1;
    bool is_tracking = false;
    double time = 0.0;
    std::vector<double> x = {state.x};
    std::vector<double> y = {state.y};
    std::vector<double> yaw = {state.yaw};
    std::vector<double> v = {state.v};
    std::vector<double> t = {0.0};

    std::cout << "Initial target index: " << target_idx << ", coordinates: (" << cx[target_idx] << ", " << cy[target_idx] << ")" << std::endl;

    while (max_simulation_time >= time && (!is_tracking || last_idx > target_idx)) {
        double ai = pid_control(target_speed, state.v);

        if (!is_tracking) {
            target_idx = calc_target_index(state, cx, cy).first;
            double delta = simple_steering_control(state, cx[target_idx], cy[target_idx]);
            state.update(ai, delta);

            double dist_to_target = state.distance_to(cx[target_idx], cy[target_idx]);
            if (dist_to_target < approach_threshold) {
                is_tracking = true;
                std::cout << "Reached target " << target_idx << " (" << cx[target_idx] << ", " << cy[target_idx] 
                          << "), starting path tracking" << std::endl;
            }
        } else {
            auto [di, new_target_idx] = stanley_control(state, cx, cy, cyaw, target_idx);
            target_idx = new_target_idx;
            state.update(ai, di);
        }

        time += dt;
        x.push_back(state.x);
        y.push_back(state.y);
        yaw.push_back(state.yaw);
        v.push_back(state.v);
        t.push_back(time);

        if (show_animation) {
            plt::cla();
            std::map<std::string, std::string> keywords_course = {{"label", "course"}, {"color", "red"}, {"marker", "."}};
            plt::plot(cx, cy, keywords_course);
            std::map<std::string, std::string> keywords_traj = {{"label", "trajectory"}, {"color", "blue"}, {"linestyle", "-"}};
            plt::plot(x, y, keywords_traj);
            std::map<std::string, std::string> keywords_target = {{"label", "target"}, {"marker", "x"}, {"color", "green"}};
            plt::plot({cx[target_idx]}, {cy[target_idx]}, keywords_target);
            plt::axis("equal");
            plt::grid(true);
            plt::xlabel("X [m]");
            plt::ylabel("Y [m]");
            plt::title("Speed [km/h]: " + std::to_string(state.v * 3.6).substr(0, 4) + 
                       (is_tracking ? " (Tracking)" : " (Approaching)"));
            plt::pause(0.001);
        }
    }

    if (last_idx < target_idx) {
        std::cerr << "Failed to reach target" << std::endl;
    }

    if (show_animation) {
        plt::cla();
        std::map<std::string, std::string> keywords_course = {{"label", "course"}, {"color", "red"}, {"marker", "."}};
        plt::plot(cx, cy, keywords_course);
        std::map<std::string, std::string> keywords_traj = {{"label", "trajectory"}, {"color", "blue"}, {"linestyle", "-"}};
        plt::plot(x, y, keywords_traj);
        plt::legend();
        plt::xlabel("X [m]");
        plt::ylabel("Y [m]");
        plt::axis("equal");
        plt::grid(true);

        plt::figure();
        std::map<std::string, std::string> keywords_speed = {{"color", "red"}, {"linestyle", "-"}};
        std::vector<double> speed_kmh(v.size());
        for (size_t i = 0; i < v.size(); ++i) {
            speed_kmh[i] = v[i] * 3.6;
        }
        plt::plot(t, speed_kmh, keywords_speed);
        plt::xlabel("Time [s]");
        plt::ylabel("Speed [km/h]");
        plt::grid(true);
        plt::show();
    }
}

int main() {
    std::cout << "Stanley path tracking simulation started" << std::endl;
    main_simulation();
    return 0;
}