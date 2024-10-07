// Dynamic Symbolic Expression Network (DSEN)
#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <limits>
#include <string>
#include <random>
#include <algorithm>  // For avoiding specific operations

#define RESET "\033[0m"
#define RED "\033[31m"
#define GREEN "\033[32m"
#define YELLOW "\033[33m"
#define CYAN "\033[36m"
#define BOLD "\033[1m"

// Random number generator for initial parameter guesses
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0.5, 2.0);

// Define basic math operations with parameters
double power_x_c(double x, double c) { return std::pow(x, c); }
double power_c_x(double x, double c) { return std::pow(c, x); }
double multiply_w(double x, double w) { return x * w; }
double add_b(double x, double b) { return x + b; }
double subtract_b(double x, double b) { return x - b; }
double divide_by_b(double x, double b) { return (b != 0) ? x / b : std::numeric_limits<double>::max(); }
double log_base_b(double x, double b) { return (x > 0 && b > 0 && b != 1) ? std::log(x) / std::log(b) : 0.0; }
double exp_b(double x, double b) { return std::exp(b * x); }
double sine(double x, double c) { return std::sin(c * x); }
double cosine(double x, double c) { return std::cos(c * x); }
double tangent(double x, double c) { return std::tan(c * x); }
double inverse(double x, double c) { return (x != 0) ? 1.0 / x : std::numeric_limits<double>::max(); }
double sqrt_x(double x, double c) { return std::sqrt(x); }
double cbrt_x(double x, double c) { return std::cbrt(x); }

// Define parametric operations with string labels
std::vector<std::pair<std::function<double(double, double)>, std::string>> parametric_operations = {
    {power_x_c, "x^c"}, {power_c_x, "c^x"}, {multiply_w, "*w"}, {add_b, "+b"}, {subtract_b, "-b"},
    {divide_by_b, "/b"}, {log_base_b, "log_b(x)"}, {exp_b, "exp(b * x)"}, {sine, "sin(c * x)"},
    {cosine, "cos(c * x)"}, {tangent, "tan(c * x)"}, {inverse, "1/x"}, {sqrt_x, "sqrt(x)"},
    {cbrt_x, "cbrt(x)"}
};

// Mean squared error cost function
double cost_function(const std::vector<double>& predictions, const std::vector<double>& targets) {
    double mse = 0.0;
    for (size_t i = 0; i < predictions.size(); ++i) {
        mse += std::pow(predictions[i] - targets[i], 2);
    }
    return mse / predictions.size();
}

// Function to optimize parameters for operations
double optimize_parameter(std::function<double(double, double)> op, const std::vector<double>& inputs, const std::vector<double>& targets) {
    double best_param = dis(gen);
    double best_cost = std::numeric_limits<double>::max();

    for (int i = 0; i < 100; ++i) {  // 100 iterations to optimize the parameter
        double param = dis(gen);
        std::vector<double> predictions;
        for (const auto& input : inputs) {
            predictions.push_back(op(input, param));
        }
        double cost = cost_function(predictions, targets);
        if (cost < best_cost) {
            best_cost = cost;
            best_param = param;
        }
    }
    return best_param;
}

// Select the best parametric operation while avoiding certain operators
std::pair<std::function<double(double)>, std::string> select_best_parametric_operation(
    const std::vector<double>& inputs, const std::vector<double>& targets,
    const std::vector<std::string>& avoid_operators) {

    double best_cost = std::numeric_limits<double>::max();
    std::function<double(double)> best_operation;
    std::string best_op_name;

    for (const auto& op_pair : parametric_operations) {
        // Check if the current operation is in the avoid list
        if (std::find(avoid_operators.begin(), avoid_operators.end(), op_pair.second) != avoid_operators.end()) {
            continue;
        }

        double param = optimize_parameter(op_pair.first, inputs, targets);
        std::vector<double> predictions;
        for (const auto& input : inputs) {
            predictions.push_back(op_pair.first(input, param));
        }

        double cost = cost_function(predictions, targets);
        if (cost < best_cost) {
            best_cost = cost;
            best_operation = [op = op_pair.first, param](double x) { return op(x, param); };  // Capture the optimized parameter
            best_op_name = op_pair.second + "(param=" + std::to_string(param) + ")";
        }
    }

    return { best_operation, best_op_name };
}

#include <iomanip>  // For setting precision in output

// Function to optimize parameters for operations
double optimize_parameter(std::function<double(double, double)> op, const std::vector<double>& inputs, const std::vector<double>& targets, double& best_param) {
    best_param = dis(gen); // Pass by reference to capture the best parameter
    double best_cost = std::numeric_limits<double>::max();

    for (int i = 0; i < 100; ++i) {  // 100 iterations to optimize the parameter
        double param = dis(gen);
        std::vector<double> predictions;
        for (const auto& input : inputs) {
            predictions.push_back(op(input, param));
        }
        double cost = cost_function(predictions, targets);
        if (cost < best_cost) {
            best_cost = cost;
            best_param = param;
        }
    }
    return best_cost;
}

// Select the best parametric operation while avoiding certain operators
std::pair<std::function<double(double)>, std::string> select_best_parametric_operation(
    const std::vector<double>& inputs, const std::vector<double>& targets,
    const std::vector<std::string>& avoid_operators, std::vector<double>& cost_per_op) {

    double best_cost = std::numeric_limits<double>::max();
    std::function<double(double)> best_operation;
    std::string best_op_name;

    int op_index = 0;  // To track the index of the operation
    for (const auto& op_pair : parametric_operations) {
        // Check if the current operation is in the avoid list
        if (std::find(avoid_operators.begin(), avoid_operators.end(), op_pair.second) != avoid_operators.end()) {
            cost_per_op.push_back(std::numeric_limits<double>::max());  // Mark avoided ops with high cost
            continue;
        }

        double param;
        double cost = optimize_parameter(op_pair.first, inputs, targets, param);
        cost_per_op.push_back(cost);  // Store the cost for this operation

        if (cost < best_cost) {
            best_cost = cost;
            best_operation = [op = op_pair.first, param](double x) { return op(x, param); };  // Capture the optimized parameter
            best_op_name = op_pair.second + "(param=" + std::to_string(param) + ")";
        }
        op_index++;
    }

    return { best_operation, best_op_name };
}

// Multi-layer model with parametric operations and avoid list
void multi_layer_model(std::vector<double> inputs, const std::vector<double>& targets,
    double gamma, int max_layers, double threshold, double tolerance,
    const std::vector<std::string>& avoid_operators) {

    std::vector<std::string> operations_used;
    double current_cost = cost_function(inputs, targets);
    double previous_cost = current_cost;

    // Cost matrix to store costs of each operation at each layer
    std::vector<std::vector<double>> cost_matrix;

    for (int layer = 0; layer < max_layers; ++layer) {
        std::cout << BOLD << CYAN << "\nLayer " << layer + 1 << ":" << RESET << "\n";

        std::vector<double> cost_per_op;
        // Select the best parametric operation, excluding the ones in the avoid list
        auto best_parametric_op = select_best_parametric_operation(inputs, targets, avoid_operators, cost_per_op);
        cost_matrix.push_back(cost_per_op);  // Store the costs for this layer

        std::vector<double> parametric_predictions = inputs;
        for (size_t i = 0; i < inputs.size(); ++i) {
            parametric_predictions[i] = best_parametric_op.first(inputs[i]);
        }
        double parametric_cost = cost_function(parametric_predictions, targets);

        if (parametric_cost < current_cost) {
            std::cout << GREEN << "Selected parametric operation: " << best_parametric_op.second << RESET << "\n";
            inputs = parametric_predictions;
            current_cost = parametric_cost;
            operations_used.push_back(best_parametric_op.second);
        }
        else {
            std::cout << YELLOW << "No significant improvement at layer " << layer + 1 << ". Continuing with tolerance." << RESET << "\n";
            if (std::abs(previous_cost - current_cost) < tolerance) {
                std::cout << RED << "No further improvement. Stopping..." << RESET << "\n";
                break;
            }
        }

        previous_cost = current_cost;

        std::cout << CYAN << "Cost after layer " << layer + 1 << ": " << current_cost << RESET << "\n";
        if (current_cost < threshold) {
            std::cout << GREEN << "Threshold reached. Stopping..." << RESET << "\n";
            break;
        }

        std::cout << "Outputs after layer " << layer + 1 << ": ";
        for (const auto& val : inputs) {
            std::cout << val << " ";
        }
        std::cout << "\n";
    }

    // Display the symbolic formula (sequence of operations)
    std::cout << BOLD << "Symbolic formula (operations used): " << RESET;
    for (const auto& op : operations_used) {
        std::cout << CYAN << op << " -> " << RESET;
    }
    std::cout << CYAN << "End" << RESET << "\n";

    // Display the cost matrix
    std::cout << BOLD << "\nCost matrix:" << RESET << "\n";
    std::cout << std::fixed << std::setprecision(6);
    for (size_t layer = 0; layer < cost_matrix.size(); ++layer) {
        std::cout << "Layer " << layer + 1 << ": ";
        for (const auto& cost : cost_matrix[layer]) {
            std::cout << cost << " ";
        }
        std::cout << "\n";
    }
}



int main() {
    std::vector<double> inputs = { 16.0, 36, 49.0, 100.0, 144.0  };
    std::vector<double> targets = { 4.0, 6.0, 7.0, 10.0, 12.0 };

    double gamma = 0.1;
    int max_layers = 10;
    double threshold = 0.0001;
    double tolerance = 0.001;

    // List of operators to avoid (for example, avoid sqrt and cbrt)
    std::vector<std::string> avoid_operators = { "sqrt(x)", "cbrt(x)", "x^c"};

    multi_layer_model(inputs, targets, gamma, max_layers, threshold, tolerance, avoid_operators);

    return 0;
}


