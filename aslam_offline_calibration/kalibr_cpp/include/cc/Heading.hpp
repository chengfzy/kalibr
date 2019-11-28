#pragma once
#include <ostream>
#include <string>
#include <type_traits>

namespace cc {

/**
 * @brief Heading type
 */
enum class HeadingType {
    Section,
    SubSection,
    Paragraph,
};

/**
 * @brief The default setting for whether to add a break line for heading type
 * @param type Heading type
 * @return True for add a break line, otherwise return false
 */
inline bool defaultBreakLine(HeadingType type) { return type != HeadingType::Paragraph; }

/**
 * @brief Heading to print or show some formatted text
 * @tparam Type     Heading type
 * @tparam SecLen   Print length for section and subsection type, and the length of paragraph will be SecLen / 1.5
 */
template <HeadingType Type = HeadingType::Section, unsigned short SecLen = 100>
class Heading {
  public:
    /**
     * @brief Constructor
     * @param text      Heading text
     * @param breakLine Flag to indict whether add a new line to show this info
     */
    explicit Heading(std::string&& text, bool breakLine = defaultBreakLine(Type));

    /**
     * @brief Constructor
     * @param text      Heading text
     * @param breakLine Flag to indict whether add a new line to show this info
     */
    explicit Heading(const std::string& text, bool breakLine = defaultBreakLine(Type));

    /**
     * @brief Destructor
     */
    ~Heading() = default;

  public:
    /**
     * @brief Print heading
     * @param os    Output stream
     * @param info  Heading info
     * @return Output stream
     */
    template <HeadingType Type_, unsigned short ParLen_>
    friend std::ostream& operator<<(std::ostream& os, const Heading<Type_, ParLen_>& info);

  private:
    std::string text_;  // heading text
    bool breakLine_;    // flag to indict whether add a break line to show this heading
};

/**************************************** Type Definition ****************************************/
using Section = Heading<HeadingType::Section>;
using SubSection = Heading<HeadingType::SubSection>;
using Paragraph = Heading<HeadingType::Paragraph>;
}  // namespace cc

#include "implementation/Heading.hpp"