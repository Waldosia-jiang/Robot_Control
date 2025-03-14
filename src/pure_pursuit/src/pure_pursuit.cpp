#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include <matplotlibcpp.h>

namespace plt = matplotlibcpp;

// 参数定义
const double k = 0.1;    // 前视增益
const double Lfc = 2.0;  // 前视距离 [m]
const double Kp = 1.0;   // 速度比例增益
const double dt = 0.1;   // 时间步长 [s]
const double WB = 2.9;   // 车辆轴距 [m]
const bool show_animation = true; // 是否显示动画

// 车辆状态类
class State {
public:
    double x, y, yaw, v, rear_x, rear_y; // 位置x, y，偏航角，速度，后轮位置x, y

    // 构造函数，初始化状态
    State(double x_ = 0.0, double y_ = 0.0, double yaw_ = 0.0, double v_ = 0.0)
        : x(x_), y(y_), yaw(yaw_), v(v_) {
        updateRear();
    }

    // 更新车辆状态
    void update(double a, double delta) {
        x += v * std::cos(yaw) * dt;       // 更新x坐标
        y += v * std::sin(yaw) * dt;       // 更新y坐标
        yaw += v / WB * std::tan(delta) * dt; // 更新偏航角
        v += a * dt;                       // 更新速度
        updateRear();                      // 更新后轮位置
    }

    // 计算与某点的距离
    double calcDistance(double point_x, double point_y) const {
        double dx = rear_x - point_x;
        double dy = rear_y - point_y;
        return std::hypot(dx, dy);
    }

private:
    // 更新后轮位置
    void updateRear() {
        rear_x = x - (WB / 2) * std::cos(yaw);
        rear_y = y - (WB / 2) * std::sin(yaw);
    }
};

// 存储车辆状态历史
class States {
public:
    std::vector<double> x, y, yaw, v, t; // x, y，偏航角，速度，时间

    // 添加状态记录
    void append(double time, const State& state) {
        x.push_back(state.x);
        y.push_back(state.y);
        yaw.push_back(state.yaw);
        v.push_back(state.v);
        t.push_back(time);
    }
};

// PID速度控制（仅比例控制）
double proportionalControl(double target, double current) {
    return Kp * (target - current); // 计算加速度
}

// 目标路径类
class TargetCourse {
public:
    std::vector<double> cx, cy; // 路径的x, y坐标
    int old_nearest_point_index; // 上一次最近点索引

    // 构造函数，初始化路径
    TargetCourse(const std::vector<double>& cx_, const std::vector<double>& cy_)
        : cx(cx_), cy(cy_), old_nearest_point_index(-1) {}

    // 搜索目标点索引和前视距离
    std::pair<int, double> searchTargetIndex(const State& state) {
        int ind;
        if (old_nearest_point_index == -1) { // 第一次搜索最近点
            std::vector<double> d;
            for (size_t i = 0; i < cx.size(); ++i) {
                double dx = state.rear_x - cx[i];
                double dy = state.rear_y - cy[i];
                d.push_back(std::hypot(dx, dy));
            }
            ind = std::min_element(d.begin(), d.end()) - d.begin();
            old_nearest_point_index = ind;
        } else { // 从上一次最近点开始搜索
            ind = old_nearest_point_index;
            double distance_this_index = state.calcDistance(cx[ind], cy[ind]);
            while (true) {
                double distance_next_index = state.calcDistance(cx[ind + 1], cy[ind + 1]);
                if (distance_this_index < distance_next_index) break;
                ind = (ind + 1 < static_cast<int>(cx.size())) ? ind + 1 : ind;
                distance_this_index = distance_next_index;
            }
            old_nearest_point_index = ind;
        }

        double Lf = k * state.v + Lfc; // 更新前视距离
        while (Lf > state.calcDistance(cx[ind], cy[ind])) { // 寻找前视点
            if (ind + 1 >= static_cast<int>(cx.size())) break;
            ind++;
        }
        return {ind, Lf};
    }
};

// Pure Pursuit转向控制
std::pair<double, int> purePursuitSteerControl(const State& state, TargetCourse& trajectory, int pind) {
    std::pair<int, double> target_info = trajectory.searchTargetIndex(state);
    int ind = target_info.first;    // 目标点索引
    double Lf = target_info.second; // 前视距离

    if (pind >= ind) ind = pind; // 确保索引不回退

    // 获取目标点坐标
    double tx = (ind < static_cast<int>(trajectory.cx.size())) ? trajectory.cx[ind] : trajectory.cx.back();
    double ty = (ind < static_cast<int>(trajectory.cy.size())) ? trajectory.cy[ind] : trajectory.cy.back();

    double alpha = std::atan2(ty - state.rear_y, tx - state.rear_x) - state.yaw; // 计算角度差
    double delta = std::atan2(2.0 * WB * std::sin(alpha) / Lf, 1.0); // 计算转向角

    return {delta, ind};
}

// 绘制箭头表示车辆方向
void plotArrow(double x, double y, double yaw, double length = 1.0, double width = 0.5) {
    std::vector<double> x_vec = {x, x + length * std::cos(yaw)};
    std::vector<double> y_vec = {y, y + length * std::sin(yaw)};
    plt::plot(x_vec, y_vec, "r-"); // 绘制箭身
    plt::plot({x}, {y}, "ro");     // 绘制起点
}

// 主模拟函数
void mainSimulation() {
    std::vector<double> cx, cy; // 目标路径
    for (double i = 0; i < 50; i += 0.5) { // 生成正弦路径
        cx.push_back(i);
        cy.push_back(std::sin(i / 5.0) * i / 2.0);
    }

    double target_speed = 10.0 / 3.6; // 目标速度 [m/s]
    double T = 100.0;                 // 最大模拟时间

    State state(-0.0, -3.0, 0.0, 0.0); // 初始状态
    States states;
    states.append(0.0, state);         // 记录初始状态

    TargetCourse target_course(cx, cy); // 创建目标路径
    int target_ind = target_course.searchTargetIndex(state).first; // 初始目标点
    int lastIndex = cx.size() - 1;    // 路径终点索引
    double time = 0.0;                // 当前时间

    while (T >= time && lastIndex > target_ind) { // 模拟循环
        double ai = proportionalControl(target_speed, state.v); // 计算加速度
        std::pair<double, int> steer_info = purePursuitSteerControl(state, target_course, target_ind);
        double di = steer_info.first;   // 转向角
        target_ind = steer_info.second; // 更新目标点索引

        state.update(ai, di); // 更新车辆状态
        time += dt;           // 增加时间
        states.append(time, state); // 记录状态

        if (show_animation) { // 显示动画
            plt::cla();       // 清空画布
            plotArrow(state.x, state.y, state.yaw); // 绘制车辆方向
            std::map<std::string, std::string> keywords_course = {{"label", "course"}, {"color", "red"}, {"linestyle", "-"}};
            plt::plot(cx, cy, keywords_course); // 绘制目标路径
            std::map<std::string, std::string> keywords_traj = {{"label", "trajectory"}, {"color", "blue"}, {"linestyle", "-"}};
            plt::plot(states.x, states.y, keywords_traj); // 绘制实际轨迹
            std::map<std::string, std::string> keywords_target = {{"label", "target"}, {"marker", "x"}, {"color", "green"}};
            plt::plot({cx[target_ind]}, {cy[target_ind]}, keywords_target); // 绘制目标点
            plt::axis("equal"); // 设置坐标轴等比例
            plt::grid(true);    // 显示网格
            plt::title("速度[km/h]: " + std::to_string(state.v * 3.6).substr(0, 4)); // 显示速度
            plt::pause(0.001);  // 暂停以显示动画
        }
    }

    if (show_animation) { // 显示最终结果
        plt::cla();       // 清空画布
        std::map<std::string, std::string> keywords_course = {{"label", "course"}, {"color", "red"}, {"marker", "."}};
        plt::plot(cx, cy, keywords_course); // 绘制目标路径
        std::map<std::string, std::string> keywords_traj = {{"label", "trajectory"}, {"color", "blue"}, {"linestyle", "-"}};
        plt::plot(states.x, states.y, keywords_traj); // 绘制实际轨迹
        plt::legend();    // 显示图例
        plt::xlabel("x[m]"); // x轴标签
        plt::ylabel("y[m]"); // y轴标签
        plt::axis("equal");  // 设置坐标轴等比例
        plt::grid(true);     // 显示网格

        plt::figure();       // 新建画布
        std::map<std::string, std::string> keywords_speed = {{"color", "red"}, {"linestyle", "-"}};
        std::vector<double> speed_kmh(states.v.size());
        for (size_t i = 0; i < states.v.size(); ++i) speed_kmh[i] = states.v[i] * 3.6; // 转换为km/h
        plt::plot(states.t, speed_kmh, keywords_speed); // 绘制速度曲线
        plt::xlabel("时间[s]"); // x轴标签
        plt::ylabel("速度[km/h]"); // y轴标签
        plt::grid(true);     // 显示网格
        plt::show();         // 显示所有图形
    }
}

// 主函数
int main() {
    std::cout << "Pure pursuit路径跟踪模拟开始" << std::endl;
    mainSimulation(); // 运行模拟
    return 0;
}
