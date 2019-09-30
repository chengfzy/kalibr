#pragma once
#include <iomanip>
#include <iostream>
#include <string>
#include "Eigen/Core"

/**
 * @brief Debug information type
 */
enum DebugInfoType {
    Section,
    SubSection,
    Paragraph,
};


// print debug information
void printInfo(const std::string& title, DebugInfoType type = DebugInfoType::Section, bool breakLine = true);

/**
 * @brief Print sparse matrix
 */
void printSpareMatrix(const Eigen::MatrixXd& m, std::ostream& os);
