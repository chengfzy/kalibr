namespace cc {

// Constructor
template <HeadingType Type, unsigned short Len>
Heading<Type, Len>::Heading(std::string&& text, bool breakLine) : text_(std::move(text)), breakLine_(breakLine) {}

// Constructor
template <HeadingType Type, unsigned short Len>
Heading<Type, Len>::Heading(const std::string& text, bool breakLine) : text_(text), breakLine_(breakLine) {}

// Print heading
template <HeadingType Type, unsigned short Len>
std::ostream& operator<<(std::ostream& os, const Heading<Type, Len>& info) {
    if (info.breakLine_) {
        os << std::endl;
    }
    char fillChar('=');
    switch (Type) {
        case HeadingType::Section:
            fillChar = '=';
            break;
        case HeadingType::SubSection:
            fillChar = '*';
            break;
        case HeadingType::Paragraph:
            fillChar = '-';
            break;
    }

    if (info.text_.empty()) {
        os << std::string(Len, fillChar);
    } else {
        std::string fillStr(std::max(5, static_cast<int>((Len - info.text_.size() - 1) / 2)), fillChar);
        os << fillStr << " " << info.text_ << " " << fillStr;
    }
    os << std::endl;

    return os;
}

}  // namespace cc