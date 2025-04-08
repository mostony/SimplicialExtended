#pragma once

#include <vector>

#include "AbstractModel.h"

class CombinatorialComplex : public AbstractModel {
   public:
    void Build(std::vector<std::vector<int>> data);

    void BuildWithRank(std::vector<std::vector<int>> data, std::vector<int> ranks);

    std::vector<std::vector<int>> GetSubsets();
    void BuildFromBinary(std::vector<std::vector<int>> binary, bool on_column) override;
};
