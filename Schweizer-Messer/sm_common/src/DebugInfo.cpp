#include "sm/DebugInfo.h"

#define CC_DEBUG true

// print debug information
void printInfo(const std::string& title, DebugInfoType type, bool breakLine) {
#if defined(CC_DEBUG)
    if (breakLine) {
        std::cout << std::endl;
    }

    // info length and fill char
    int infoLen = 0;
    char fillChar = '=';
    switch (type) {
        case DebugInfoType::Section:
            infoLen = 100;
            fillChar = '=';
            break;
        case DebugInfoType::SubSection:
            infoLen = 100;
            fillChar = '*';
            break;
        case DebugInfoType::Paragraph:
            infoLen = 67;
            fillChar = '-';
            break;
    }

    if (title.empty()) {
        std::cout << std::string(infoLen, fillChar);
    } else {
        std::string fillStr(std::max(5, static_cast<int>((infoLen - title.size() - 1) / 2)), fillChar);
        std::cout << fillStr << " " << title << " " << fillStr;
    }
    std::cout << std::endl;
#endif
}

/**
 * @brief Print sparse matrix
 */
void printSpareMatrix(const Eigen::MatrixXd& m, std::ostream& os) {
    for (int i = 0; i < m.rows(); ++i) {
        for (int j = 0; j < m.cols(); ++j) {
            if (std::abs(m(i, j)) > std::numeric_limits<double>::epsilon()) {
                os << i << ", " << j << ", " << std::setprecision(10) << m(i, j) << std::endl;
            }
        }
    }
}