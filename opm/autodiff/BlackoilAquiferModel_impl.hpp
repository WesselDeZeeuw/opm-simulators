namespace Opm {

  template<typename TypeTag>
  BlackoilAquiferModel<TypeTag>::
  BlackoilAquiferModel(Simulator& ebosSimulator)
  : ebosSimulator_(ebosSimulator)
  {
    init();
  }


  // called at the end of a time step
  template<typename TypeTag>
  void
  BlackoilAquiferModel<TypeTag>:: timeStepSucceeded(const SimulatorTimerInterface& timer)
  {
    if(aquiferCarterTracyActive())

    {
      for (auto aquifer = aquifers_CarterTracy.begin(); aquifer != aquifers_CarterTracy.end(); ++aquifer)
      {
        aquifer->afterTimeStep(timer);
      }
    }
    if(aquiferFetkovichActive())
    {
      for (auto aquifer = aquifers_Fetkovich.begin(); aquifer != aquifers_Fetkovich.end(); ++aquifer)
      {
        aquifer->afterTimeStep(timer);
      }

    }
  }

  template<typename TypeTag>
  void
  BlackoilAquiferModel<TypeTag>::
  assemble( const SimulatorTimerInterface& timer,
    const int iterationIdx                )
    {
      if ( !aquiferActive() ) {
        return;
      }

      // We need to update the reservoir pressures connected to the aquifer
      updateConnectionIntensiveQuantities();

      if (iterationIdx == 0) {
        // We can do the Table check and coefficients update in this function
        // For now, it does nothing!
        prepareTimeStep(timer);
      }

      assembleAquiferEq(timer);
    }

    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>::beginIteration()
    { }

    template<typename TypeTag>
    template <class Context>
    void
    BlackoilAquiferModel<TypeTag>::addToSource(RateVector& rates,
                                               const Context& context,
                                               unsigned spaceIdx,
                                               unsigned timeIdx) const
    {
      ElementContext elemCtx(ebosSimulator_);
      const auto& gridView = ebosSimulator_.gridView();
      const auto& elemEndIt = gridView.template end</*codim=*/0, Dune::Interior_Partition>();
      for (auto elemIt = gridView.template begin</*codim=*/0, Dune::Interior_Partition>();
      elemIt != elemEndIt;
      ++elemIt)
      {
        elemCtx.updatePrimaryStencil(*elemIt);
        elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
      }
    }

    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: assembleAquiferEq(const SimulatorTimerInterface& timer)
    {
      if(aquiferCarterTracyActive()){
        for (auto aquifer = aquifers_CarterTracy.begin(); aquifer != aquifers_CarterTracy.end(); ++aquifer)
        {
          aquifer->assembleAquiferEq(timer);
        }
      }
      if(aquiferFetkovichActive()){
        for (auto aquifer = aquifers_Fetkovich.begin(); aquifer != aquifers_Fetkovich.end(); ++aquifer)
        {
          aquifer->assembleAquiferEq(timer);
        }
      }
      else if (!aquiferCarterTracyActive())
      {
        for (auto aquifer = aquifers_Fetkovich.begin(); aquifer != aquifers_Fetkovich.end(); ++aquifer)
        {
          aquifer->assembleAquiferEq(timer);
        }
      }
    }

    template<typename TypeTag>
    void BlackoilAquiferModel<TypeTag>:: prepareTimeStep(const SimulatorTimerInterface& timer)
    {
      if(aquiferCarterTracyActive())
      {
        for (auto aquifer = aquifers_CarterTracy.begin(); aquifer != aquifers_CarterTracy.end(); ++aquifer)
        {
          aquifer->beforeTimeStep(timer);
        }
      }
      if(aquiferFetkovichActive())
      {
        for (auto aquifer = aquifers_Fetkovich.begin(); aquifer != aquifers_Fetkovich.end(); ++aquifer)
        {
          aquifer->beforeTimeStep(timer);
        }
      }
      else if (!aquiferCarterTracyActive())
      {
        for (auto aquifer = aquifers_Fetkovich.begin(); aquifer != aquifers_Fetkovich.end(); ++aquifer)
        {
          aquifer->beforeTimeStep(timer);
        }
      }
    }

    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>::endEpisode()
    { }

    // Initialize the aquifers in the deck
    template<typename TypeTag>
    void
    BlackoilAquiferModel<TypeTag>:: init()
    {
      const auto& deck = ebosSimulator_.vanguard().deck();
      if (deck.hasKeyword("AQUCT")) {

        updateConnectionIntensiveQuantities();
        const auto& eclState = ebosSimulator_.vanguard().eclState();

        // Get all the carter tracy aquifer properties data and put it in aquifers vector
        const AquiferCT aquiferct = AquiferCT(eclState,deck);
        const Aquancon aquifer_connect = Aquancon(eclState.getInputGrid(), deck);

        std::vector<AquiferCT::AQUCT_data> aquifersData = aquiferct.getAquifers();
        std::vector<Aquancon::AquanconOutput> aquifer_connection = aquifer_connect.getAquOutput();

        assert( aquifersData.size() == aquifer_connection.size() );


        for (size_t i = 0; i < aquifersData.size(); ++i)
        {
          aquifers_CarterTracy.push_back(
            AquiferCarterTracy<TypeTag> (aquifersData.at(i), aquifer_connection.at(i), ebosSimulator_)
          );
        }
      }
      if(deck.hasKeyword("AQUFETP"))
      {

        updateConnectionIntensiveQuantities();
        const auto& eclState = ebosSimulator_.vanguard().eclState();

        // Get all the carter tracy aquifer properties data and put it in aquifers vector
        const Aquifetp aquifetp = Aquifetp(deck);
        const Aquancon aquifer_connect = Aquancon(eclState.getInputGrid(), deck);

        std::vector<Aquifetp::AQUFETP_data> aquifersData = aquifetp.getAquifers();
        std::vector<Aquancon::AquanconOutput> aquifer_connection = aquifer_connect.getAquOutput();

        assert( aquifersData.size() == aquifer_connection.size() );
        for (size_t i = 0; i < aquifersData.size(); ++i)
        {
          aquifers_Fetkovich.push_back(
            AquiferFetkovich<TypeTag> (aquifersData.at(i), aquifer_connection.at(i), ebosSimulator_)
          );
        }
      }
    }

    template<typename TypeTag>
    bool
    BlackoilAquiferModel<TypeTag>:: aquiferActive() const
    {
      return (aquiferCarterTracyActive() || aquiferFetkovichActive());
    }
    template<typename TypeTag>
    bool
    BlackoilAquiferModel<TypeTag>:: aquiferCarterTracyActive() const
    {
      return !aquifers_CarterTracy.empty();
    }
    template<typename TypeTag>
    bool
    BlackoilAquiferModel<TypeTag>:: aquiferFetkovichActive() const
    {
      return !aquifers_Fetkovich.empty();
    }

  } // namespace Opm
