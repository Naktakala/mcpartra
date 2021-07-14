#ifndef MCPARTRA_MCPT_CUSTOMVOLUMETALLY_H
#define MCPARTRA_MCPT_CUSTOMVOLUMETALLY_H

namespace mcpartra
{

//#########################################################
/**Custom tally Structure.*/
  class CustomVolumeTally
  {
  public:
    struct TallyFluctuationChart
    {
      std::vector<double> average;
      std::vector<double> sigma;
    };
    typedef TallyFluctuationChart TFC;
  public:
    GridTallyBlock        grid_tally;
    std::vector<bool>     local_cell_tally_mask; ///< Indicates whether a local cell is part of tally
    double                tally_volume=0.0;
    std::vector<TFC>      tally_fluctation_chart;

    explicit
    CustomVolumeTally(std::vector<bool>& in_masking) :
      local_cell_tally_mask(in_masking) {}

    void Initialize(size_t tally_size, double in_volume)
    {
      grid_tally.Resize(tally_size);
      tally_volume = in_volume;
    }
  };
}//namespace mcpartra
#endif //MCPARTRA_MCPT_CUSTOMVOLUMETALLY_H
