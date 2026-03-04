#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>
#include <chrono>

using namespace std;

// Constants
const double G = 9.81;  // gravity (m/s²)
const double DT = 0.01; // time step (s)

struct State {
    double theta;      // angle (radians)
    double theta_dot;  // angular velocity (rad/s)
};

class PendulumSimulator {
private:
    double L;           // length (m)
    double damping;     // damping coefficient
    double noise_std;   // measurement noise standard deviation
    
    // RNG for noise
    mt19937 rng;
    normal_distribution<double> noise_distribution;
    
public:
    PendulumSimulator(double length = 1.0, double damp = 0.0, double noise = 0.0)
        : L(length), damping(damp), noise_std(noise), 
          rng(chrono::steady_clock::now().time_since_epoch().count()),
          noise_distribution(0.0, noise) {}
    
    // RK4 integration step
    State rk4_step(const State& state, double dt) {
        auto derivatives = [this](const State& s) -> State {
            double theta_ddot = -(G / L) * sin(s.theta) - damping * s.theta_dot;
            return {s.theta_dot, theta_ddot};
        };
        
        State k1 = derivatives(state);
        
        State temp2 = {state.theta + 0.5 * dt * k1.theta, 
                       state.theta_dot + 0.5 * dt * k1.theta_dot};
        State k2 = derivatives(temp2);
        
        State temp3 = {state.theta + 0.5 * dt * k2.theta,
                       state.theta_dot + 0.5 * dt * k2.theta_dot};
        State k3 = derivatives(temp3);
        
        State temp4 = {state.theta + dt * k3.theta,
                       state.theta_dot + dt * k3.theta_dot};
        State k4 = derivatives(temp4);
        
        State next;
        next.theta = state.theta + (dt / 6.0) * (k1.theta + 2*k2.theta + 2*k3.theta + k4.theta);
        next.theta_dot = state.theta_dot + (dt / 6.0) * (k1.theta_dot + 2*k2.theta_dot + 2*k3.theta_dot + k4.theta_dot);
        
        return next;
    }
    
    // Add measurement noise to state
    State add_noise(const State& state) {
        if (noise_std == 0.0) return state;
        
        State noisy;
        noisy.theta = state.theta + noise_distribution(rng);
        noisy.theta_dot = state.theta_dot + noise_distribution(rng);
        return noisy;
    }
    
    // Simulate single trajectory
    vector<pair<State, State>> simulate_trajectory(
        double theta0, 
        double theta_dot0,
        double duration,
        int num_steps
    ) {
        vector<pair<State, State>> trajectory;
        State current = {theta0, theta_dot0};
        
        for (int i = 0; i < num_steps - 1; i++) {
            State with_noise = add_noise(current);
            State next = rk4_step(current, DT);
            State next_with_noise = add_noise(next);
            
            trajectory.push_back({with_noise, next_with_noise});
            current = next;
        }
        
        return trajectory;
    }
    
    // Generate large dataset
    vector<pair<State, State>> generate_dataset(
        int num_trajectories,
        int steps_per_trajectory,
        double angle_min = -M_PI/4,
        double angle_max = M_PI/4,
        double vel_min = -1.0,
        double vel_max = 1.0
    ) {
        vector<pair<State, State>> all_data;
        
        uniform_real_distribution<double> angle_dist(angle_min, angle_max);
        uniform_real_distribution<double> vel_dist(vel_min, vel_max);
        
        cout << "Generating " << num_trajectories << " trajectories..." << endl;
        
        for (int i = 0; i < num_trajectories; i++) {
            if ((i + 1) % 100 == 0) {
                cout << "  Completed " << (i + 1) << " trajectories" << endl;
            }
            
            double theta0 = angle_dist(rng);
            double theta_dot0 = vel_dist(rng);
            
            auto traj = simulate_trajectory(theta0, theta_dot0, 10.0, steps_per_trajectory);
            all_data.insert(all_data.end(), traj.begin(), traj.end());
        }
        
        cout << "Total samples generated: " << all_data.size() << endl;
        return all_data;
    }
    
    // Save to CSV
    void save_to_csv(const vector<pair<State, State>>& data, const string& filename) {
        ofstream file(filename);
        file << "angle_t,velocity_t,angle_t1,velocity_t1\n";
        
        for (const auto& [current, next] : data) {
            file << fixed << setprecision(6)
                 << current.theta << ","
                 << current.theta_dot << ","
                 << next.theta << ","
                 << next.theta_dot << "\n";
        }
        
        file.close();
        cout << "Data saved to " << filename << endl;
    }
};

int main(int argc, char* argv[]) {
    // Configuration
    int num_trajectories = 1000;
    int steps_per_trajectory = 100;
    double measurement_noise = 0.05;  // 0.05 rad noise
    
    // Parse command line arguments
    if (argc > 1) num_trajectories = stoi(argv[1]);
    if (argc > 2) steps_per_trajectory = stoi(argv[2]);
    if (argc > 3) measurement_noise = stod(argv[3]);
    
    cout << "=== Pendulum Data Simulator ===" << endl;
    cout << "Trajectories: " << num_trajectories << endl;
    cout << "Steps per trajectory: " << steps_per_trajectory << endl;
    cout << "Measurement noise std: " << measurement_noise << " rad" << endl;
    cout << endl;
    
    auto start = chrono::high_resolution_clock::now();
    
    // Create simulator with small damping and noise
    PendulumSimulator sim(1.0, 0.1, measurement_noise);
    
    // Generate dataset
    auto data = sim.generate_dataset(
        num_trajectories,
        steps_per_trajectory,
        -M_PI/6, M_PI/6,  // angle range
        -2.0, 2.0          // velocity range
    );
    
    // Save to CSV
    sim.save_to_csv(data, "pendulum_data.csv");
    
    auto end = chrono::high_resolution_clock::now();
    double elapsed = chrono::duration<double>(end - start).count();
    
    cout << "\nGeneration time: " << fixed << setprecision(2) << elapsed << " seconds" << endl;
    cout << "Samples per second: " << (data.size() / elapsed) << endl;
    
    return 0;
}