#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "DataFormats/MuonData/interface/MuonDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"


// Do not work for MuonDigiCollection, because it is not derived from std::vector
//typedef SingleObjectSelector<
//          CSCCorrelatedLCTDigiCollection,
//          StringCutObjectSelector<CSCCorrelatedLCTDigi>
//        > MuonLCTDigiSelector;

// Have to make my own version
// Adapted from CommonTools/UtilAlgos/interface/SelectedOutputCollectionTrait.h
namespace helper {
  template<typename IndexType, typename DigiType>
  struct SelectedOutputCollectionTrait<MuonDigiCollection<IndexType, DigiType> > {
    typedef typename std::vector<DigiType> type;
  };
}

// Adapted from CommonTools/UtilAlgos/interface/SingleElementCollectionSelector.h
template<typename InputCollection, typename Selector,
     typename OutputCollection = typename ::helper::SelectedOutputCollectionTrait<InputCollection>::type,
     typename StoreContainer = typename ::helper::StoreContainerTrait<OutputCollection>::type,
     typename RefAdder = typename ::helper::SelectionAdderTrait<InputCollection, StoreContainer>::type>
struct SingleElementMuonDigiCollectionSelector {
  typedef InputCollection collection;
  typedef StoreContainer container;
  typedef Selector selector;
  typedef typename container::const_iterator const_iterator;
  SingleElementMuonDigiCollectionSelector(const edm::ParameterSet & cfg, edm::ConsumesCollector && iC) :
    select_(reco::modules::make<Selector>(cfg, iC)) { }
  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }
  void select(const edm::Handle<InputCollection> & c, const edm::Event &, const edm::EventSetup &) {
    selected_.clear();
    // Loop over a std::map<index, std::vector<digi> >
    size_t idx = 0;
    auto chamber = c->begin();
    auto chend   = c->end();
    for (; chamber != chend; ++chamber) {
      auto digi = (*chamber).second.first;
      auto dend = (*chamber).second.second;
      for (; digi != dend; ++digi) {
        if (select_(*digi)) {
          selected_.push_back(&(*digi));
        }
        ++idx;
      }
    }
  }
private:
  container selected_;
  selector select_;
  RefAdder addRef_;
  friend struct reco::modules::SingleElementCollectionSelectorEventSetupInit<SingleElementMuonDigiCollectionSelector>;
};

// Adapted from CommonTools/UtilAlgos/interface/SingleObjectSelector.h
template<typename InputCollection, typename Selector,
    typename EdmFilter = edm::EDFilter,
    typename OutputCollection = typename ::helper::SelectedOutputCollectionTrait<InputCollection>::type,
    typename StoreContainer = typename ::helper::StoreContainerTrait<OutputCollection>::type,
    typename PostProcessor = ::helper::NullPostProcessor<OutputCollection, EdmFilter>,
    typename StoreManager = typename ::helper::StoreManagerTrait<OutputCollection, EdmFilter>::type,
    typename Base = typename ::helper::StoreManagerTrait<OutputCollection, EdmFilter>::base,
    typename RefAdder = typename ::helper::SelectionAdderTrait<InputCollection, StoreContainer>::type>
class SingleMuonDigiSelector :
  public ObjectSelector<SingleElementMuonDigiCollectionSelector<InputCollection, Selector, OutputCollection, StoreContainer, RefAdder>,
      OutputCollection, NonNullNumberSelector, PostProcessor, StoreManager, Base> {
public:
  explicit SingleMuonDigiSelector( const edm::ParameterSet & cfg ) :
    ObjectSelector<SingleElementMuonDigiCollectionSelector<InputCollection, Selector, OutputCollection, StoreContainer, RefAdder>,
    OutputCollection, NonNullNumberSelector, PostProcessor, StoreManager, Base>( cfg ) { }
  virtual ~SingleMuonDigiSelector() { }
};

typedef SingleMuonDigiSelector<
          CSCCorrelatedLCTDigiCollection,
          StringCutObjectSelector<CSCCorrelatedLCTDigi>
        > MuonLCTDigiSelector;

DEFINE_FWK_MODULE(MuonLCTDigiSelector);
