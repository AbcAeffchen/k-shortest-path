#ifndef SRC_TOOLS_OUTPUT
#define SRC_TOOLS_OUTPUT

#include <iostream>
#include <concepts>
#include <ranges>
#include <vector>
#include "KspBasics.h"

class Print
{
public:
    static auto& data() noexcept
    {
        return std::cout;
    }

    static auto& info() noexcept
    {
        return std::cerr;
    }
};

class JSON
{
private:
    std::string jsonStr;

    void newItem() noexcept
    {
        if(!jsonStr.empty())
            jsonStr += ",";
    }

    void writeName(const std::string& name) noexcept
    {
        jsonStr += "\"" + name + "\":";
    }

    template<typename T>
    requires ((std::is_floating_point_v<T> || std::is_integral_v<T>) && (!std::is_same_v<T, bool>))
    void writeData(const T& data) noexcept
    {
        jsonStr += std::to_string(data);
    }

    void writeData(const bool data) noexcept
    {
        if(data)
            jsonStr += "true";
        else
            jsonStr += "false";
    }

    void writeData(const std::string& data) noexcept
    {
        assert(data.find('\"') == std::string::npos);

        jsonStr += "\"" + data + "\"";
    }

    void writeData(const JSON& data) noexcept
    {
        jsonStr += data.getString();
    }

    template<typename WeightType>
    void writeData(const KSPPath<WeightType>& data) noexcept
    {
        JSON json;
        json.add("path", data.path);
        json.add("length", data.length);
        json.add("deviationNodeIndex", data.deviationNodeIndex);
        json.add("parentPathId", data.parentPathId);

        writeData(json);
    }

    template<typename T>
    void writeData(const std::vector<T>& dataList) noexcept
    {
        jsonStr += "[";
        writeData(dataList.front());

        for(const auto& item: dataList | std::views::drop(1))
        {
            jsonStr += ",";
            writeData(item);
        }

        jsonStr += "]";
    }

public:
    template<typename T>
    void add(const std::string& name, const T& data) noexcept
    {
        newItem();
        writeName(name);
        writeData(data);
    }

    [[nodiscard]] std::string getString() const noexcept
    {
        return "{" + jsonStr + "}";
    }
};

inline std::ostream& operator<<(std::ostream& os, const JSON& json) noexcept
{
    os << json.getString();
    return os;
}

#endif //SRC_TOOLS_OUTPUT
