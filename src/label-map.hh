
#ifndef LABEL_MAP_HH
#define LABEL_MAP_HH

#include <map>
#include <vector>
#include <stdexcept>
#include <string>

class LabelMap {
    size_t 							            _count;
    std::map<std::string, size_t> 	_map;
    std::vector<std::string>     	  _names;

public:
    LabelMap() : _count(0) {};

    void clear(){_count = 0; _map.clear(); _names.clear();};

    struct AlreadyPushedEx {
  		std::string label;
  		AlreadyPushedEx(std::string l) : label(l) {}
    };
    
    struct UnkownLabelEx {
  		std::string label;
  		UnkownLabelEx(std::string l) : label(l) {}
    };

    size_t push(std::string label) throw(AlreadyPushedEx);
    size_t size() const { return _count; }

    size_t operator[](std::string label) const throw(UnkownLabelEx);
    std::string name(unsigned int idx)   const throw(std::out_of_range);
};

#endif // LABEL_MAP_HH
