#include "output_connector_data_transfer/output_connection.h"

template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename FunctionSpaceType2, int nComponents2a, int nComponents2b>
void OutputConnection::
initialize(const Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b> &transferableSolutionData1,
           const Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b> &transferableSolutionData2)
{
  if (transferDirectionTerm1To2_)
  {
    nFieldVariablesTerm1Vector1_ = transferableSolutionData1.variable1.size();
    nFieldVariablesTerm1Vector2_ = transferableSolutionData1.variable2.size();
    nFieldVariablesTerm2Vector1_ = transferableSolutionData2.variable1.size();
    nFieldVariablesTerm2Vector2_ = transferableSolutionData2.variable2.size();
  }
  else
  {
    nFieldVariablesTerm1Vector1_ = transferableSolutionData2.variable1.size();
    nFieldVariablesTerm1Vector2_ = transferableSolutionData2.variable2.size();
    nFieldVariablesTerm2Vector1_ = transferableSolutionData1.variable1.size();
    nFieldVariablesTerm2Vector2_ = transferableSolutionData1.variable2.size();
  }

  if (!fieldVariableNamesInitialized_)
  {
    fieldVariableNamesInitialized_ = true;

    // collect field variable names for debugging
    for (const Data::ComponentOfFieldVariable<FunctionSpaceType1,nComponents1a> &entry : transferableSolutionData1.variable1)
    {
      std::stringstream name;
      if (entry.values)
      {
        name << entry.values->name() << "." << entry.componentNo;
      }
      else
      {
        name << "(null)." << entry.componentNo;
      }
      fieldVariableNamesTerm1Vector1_.push_back(name.str());
    }

    for (const Data::ComponentOfFieldVariable<FunctionSpaceType1,nComponents1b> &entry : transferableSolutionData1.variable2)
    {
      std::stringstream name;
      if (entry.values)
      {
        name << entry.values->name() << "." << entry.componentNo;
      }
      else
      {
        name << "(null)." << entry.componentNo;
      }
      fieldVariableNamesTerm1Vector2_.push_back(name.str());
    }
    for (const Data::ComponentOfFieldVariable<FunctionSpaceType2,nComponents2a> &entry : transferableSolutionData2.variable1)
    {
      std::stringstream name;
      if (entry.values)
      {
        name << entry.values->name() << "." << entry.componentNo;
      }
      else
      {
        name << "(null)." << entry.componentNo;
      }
      fieldVariableNamesTerm2Vector1_.push_back(name.str());
    }

    for (const Data::ComponentOfFieldVariable<FunctionSpaceType2,nComponents2b> &entry : transferableSolutionData2.variable2)
    {
      std::stringstream name;
      if (entry.values)
      {
        name << entry.values->name() << "." << entry.componentNo;
      }
      else
      {
        name << "(null)." << entry.componentNo;
      }
      fieldVariableNamesTerm2Vector2_.push_back(name.str());
    }

    initializeSlotInformation(transferableSolutionData1, transferableSolutionData2);
  }
}

template<typename FunctionSpaceType1, int nComponents1a, int nComponents1b, typename FunctionSpaceType2, int nComponents2a, int nComponents2b>
void OutputConnection::
initializeSlotInformation(const Data::OutputConnectorData<FunctionSpaceType1,nComponents1a,nComponents1b> &transferableSolutionData1,
                          const Data::OutputConnectorData<FunctionSpaceType2,nComponents2a,nComponents2b> &transferableSolutionData2)
{
  if (slotInformationInitialized_)
    return;

  LOG(DEBUG) << "initializeSlotInformation";

  // variable1 from 1 to 2
  // ----------------------
  transferDirectionTerm1To2_ = true;

  // determine if for this variable there are copy and non-copy entries. If so, non-copy is not possible and all have to be set to copy.
  bool copyRequired = false;

  // iterate over the first vector of variables
  for (int i = 0; i < transferableSolutionData1.variable1.size(); i++)
  {
    LOG(DEBUG) << "i=" << i;

    // call getSlotInformation, which outputs the warnings
    int fromVectorNo = 0;
    int fromVectorIndex = i;
    int toVectorNo = 0;
    int toVectorIndex = 0;
    bool avoidCopyIfPossible = true;
    bool slotIsConnected = getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo, toVectorIndex, avoidCopyIfPossible);

    if (!slotIsConnected)
      continue;

    if (!avoidCopyIfPossible)
    {
      copyRequired = true;
    }

    // check if other entries of the variable map to a different variable
    if (avoidCopyIfPossible)
    {
      // loop over other connections of this variable
      for (int i2 = i+1; i2 < transferableSolutionData1.variable1.size(); i2++)
      {
        LOG(DEBUG) << "i2=" << i2;

        int fromVectorNo = 0;
        int fromVectorIndex = i2;
        int toVectorNo2 = 0;
        int toVectorIndex2 = 0;
        bool avoidCopyIfPossible = true;
        bool slotIsConnected = getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo2, toVectorIndex2, avoidCopyIfPossible);

        if (!slotIsConnected)
          continue;

        if (toVectorNo2 != toVectorNo || !avoidCopyIfPossible)
        {
          LOG(DEBUG) << "slot " << i << " which is on variable1 maps to variable " << toVectorNo+1 << ", but slot "
            << i2 << ", which is also on variable1 maps to variable " << toVectorNo2+1 << ". avoidCopyIfPossible = " << avoidCopyIfPossible;
          LOG(DEBUG) << "Therefore copyRequired.";
          copyRequired = true;
        }
      }
    }
  }

  if (copyRequired)
  {
    LOG(DEBUG) << "transferDirectionTerm1To2_=" << transferDirectionTerm1To2_
      << ", copyRequired was set, set all connections for variable1 to copy";

    // set in all entries of the first field variable avoidCopyIfPossible=false
    for (int i = 0; i < transferableSolutionData1.variable1.size(); i++)
    {
      LOG(DEBUG) << " i=" << i;

      int fromVectorNo = 0;
      int fromVectorIndex = i;
      int toVectorNo = 0;
      int toVectorIndex = 0;
      bool avoidCopyIfPossible = true;
      bool slotIsConnected = getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo, toVectorIndex, avoidCopyIfPossible);

      if (!slotIsConnected)
        continue;

      int fromIndex = fromVectorNo*nFieldVariablesTerm1Vector1_ + fromVectorIndex;
      connectorTerm1To2_[fromIndex].avoidCopyIfPossible = false;
    }
  }

  // variable2 from 1 to 2
  // ----------------------
  transferDirectionTerm1To2_ = true;

  // determine if for this variable there are copy and non-copy entries. If so, non-copy is not possible and all have to be set to copy.
  copyRequired = false;

  // iterate over the first vector of variables
  for (int i = 0; i < transferableSolutionData1.variable2.size(); i++)
  {
    LOG(DEBUG) << "i=" << i;

    // call getSlotInformation, which outputs the warnings
    int fromVectorNo = 1;
    int fromVectorIndex = i;
    int toVectorNo = 0;
    int toVectorIndex = 0;
    bool avoidCopyIfPossible = true;
    bool slotIsConnected = getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo, toVectorIndex, avoidCopyIfPossible);

    if (!slotIsConnected)
      continue;

    if (!avoidCopyIfPossible)
    {
      copyRequired = true;
    }

    // check if other entries of the variable map to a different variable
    if (avoidCopyIfPossible)
    {
      // loop over other connections of this variable
      for (int i2 = i+1; i2 < transferableSolutionData1.variable2.size(); i2++)
      {
        LOG(DEBUG) << "i2=" << i2;

        int fromVectorNo = 1;
        int fromVectorIndex = i2;
        int toVectorNo2 = 0;
        int toVectorIndex2 = 0;
        bool avoidCopyIfPossible = true;
        bool slotIsConnected = getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo2, toVectorIndex2, avoidCopyIfPossible);

        if (!slotIsConnected)
          continue;

        if (toVectorNo2 != toVectorNo || !avoidCopyIfPossible)
        {
          LOG(DEBUG) << "slot " << i << " which is on variable2 maps to variable " << toVectorNo+1 << ", but slot "
            << i2 << ", which is also on variable2 maps to variable " << toVectorNo2+1 << ". avoidCopyIfPossible = " << avoidCopyIfPossible;
          LOG(DEBUG) << "Therefore copyRequired.";
          copyRequired = true;
        }
      }
    }
  }

  if (copyRequired)
  {
    LOG(DEBUG) << "transferDirectionTerm1To2_=" << transferDirectionTerm1To2_
      << ", copyRequired was set, set all connections for variable2 to copy";

    // set in all entries of the first field variable avoidCopyIfPossible=false
    for (int i = 0; i < transferableSolutionData1.variable2.size(); i++)
    {
      LOG(DEBUG) << " i=" << i;

      int fromVectorNo = 1;
      int fromVectorIndex = i;
      int toVectorNo = 0;
      int toVectorIndex = 0;
      bool avoidCopyIfPossible = true;
      bool slotIsConnected = getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo, toVectorIndex, avoidCopyIfPossible);

      if (!slotIsConnected)
        continue;

      int fromIndex = fromVectorNo*nFieldVariablesTerm1Vector1_ + fromVectorIndex;
      connectorTerm1To2_[fromIndex].avoidCopyIfPossible = false;
    }
  }


  // variable1 from 2 to 1
  // ----------------------
  transferDirectionTerm1To2_ = false;

  // determine if for this variable there are copy and non-copy entries. If so, non-copy is not possible and all have to be set to copy.
  copyRequired = false;

  // iterate over the first vector of variables
  for (int i = 0; i < transferableSolutionData2.variable1.size(); i++)
  {
    LOG(DEBUG) << "i=" << i;

    // call getSlotInformation, which outputs the warnings
    int fromVectorNo = 0;
    int fromVectorIndex = i;
    int toVectorNo = 0;
    int toVectorIndex = 0;
    bool avoidCopyIfPossible = true;
    bool slotIsConnected = getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo, toVectorIndex, avoidCopyIfPossible);

    if (!slotIsConnected)
      continue;

    if (!avoidCopyIfPossible)
    {
      copyRequired = true;
    }

    // check if other entries of the variable map to a different variable
    if (avoidCopyIfPossible)
    {
      // loop over other connections of this variable
      for (int i2 = i+1; i2 < transferableSolutionData2.variable1.size(); i2++)
      {
        LOG(DEBUG) << "i2=" << i2;

        int fromVectorNo = 0;
        int fromVectorIndex = i2;
        int toVectorNo2 = 0;
        int toVectorIndex2 = 0;
        bool avoidCopyIfPossible = true;
        bool slotIsConnected = getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo2, toVectorIndex2, avoidCopyIfPossible);

        if (!slotIsConnected)
          continue;

        if (toVectorNo2 != toVectorNo || !avoidCopyIfPossible)
        {
          LOG(DEBUG) << "slot " << i << " which is on variable1 maps to variable " << toVectorNo+1 << ", but slot "
            << i2 << ", which is also on variable1 maps to variable " << toVectorNo2+1 << ". avoidCopyIfPossible = " << avoidCopyIfPossible;
          LOG(DEBUG) << "Therefore copyRequired.";
          copyRequired = true;
        }
      }
    }
  }

  if (copyRequired)
  {
    LOG(DEBUG) << "transferDirectionTerm1To2_=" << transferDirectionTerm1To2_
      << ", copyRequired was set, set all connections for variable1 to copy";

    // set in all entries of the first field variable avoidCopyIfPossible=false
    for (int i = 0; i < transferableSolutionData2.variable1.size(); i++)
    {
      LOG(DEBUG) << " i=" << i;

      int fromVectorNo = 0;
      int fromVectorIndex = i;
      int toVectorNo = 0;
      int toVectorIndex = 0;
      bool avoidCopyIfPossible = true;
      bool slotIsConnected = getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo, toVectorIndex, avoidCopyIfPossible);

      if (!slotIsConnected)
        continue;

      int fromIndex = fromVectorNo*nFieldVariablesTerm2Vector1_ + fromVectorIndex;
      connectorTerm2To1_[fromIndex].avoidCopyIfPossible = false;
    }
  }


  // variable2 from 2 to 1
  // ----------------------
  transferDirectionTerm1To2_ = false;

  // determine if for this variable there are copy and non-copy entries. If so, non-copy is not possible and all have to be set to copy.
  copyRequired = false;

  // iterate over the first vector of variables
  for (int i = 0; i < transferableSolutionData2.variable2.size(); i++)
  {
    LOG(DEBUG) << "i=" << i;

    // call getSlotInformation, which outputs the warnings
    int fromVectorNo = 1;
    int fromVectorIndex = i;
    int toVectorNo = 0;
    int toVectorIndex = 0;
    bool avoidCopyIfPossible = true;
    bool slotIsConnected = getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo, toVectorIndex, avoidCopyIfPossible);

    if (!slotIsConnected)
      continue;

    if (!avoidCopyIfPossible)
    {
      copyRequired = true;
    }

    // check if other entries of the variable map to a different variable
    if (avoidCopyIfPossible)
    {
      // loop over other connections of this variable
      for (int i2 = i+1; i2 < transferableSolutionData2.variable2.size(); i2++)
      {
        LOG(DEBUG) << "i2=" << i2;

        int fromVectorNo = 1;
        int fromVectorIndex = i2;
        int toVectorNo2 = 0;
        int toVectorIndex2 = 0;
        bool avoidCopyIfPossible = true;
        bool slotIsConnected = getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo2, toVectorIndex2, avoidCopyIfPossible);

        if (!slotIsConnected)
          continue;

        if (toVectorNo2 != toVectorNo || !avoidCopyIfPossible)
        {
          LOG(DEBUG) << "slot " << i << " which is on variable2 maps to variable " << toVectorNo+1 << ", but slot "
            << i2 << ", which is also on variable2 maps to variable " << toVectorNo2+1 << ". avoidCopyIfPossible = " << avoidCopyIfPossible;
          LOG(DEBUG) << "Therefore copyRequired.";
          copyRequired = true;
        }
      }
    }
  }

  if (copyRequired)
  {
    LOG(DEBUG) << "transferDirectionTerm1To2_=" << transferDirectionTerm1To2_
      << ", copyRequired was set, set all connections for variable1 to copy";

    // set in all entries of the first field variable avoidCopyIfPossible=false
    for (int i = 0; i < transferableSolutionData2.variable2.size(); i++)
    {
      LOG(DEBUG) << " i=" << i;

      int fromVectorNo = 1;
      int fromVectorIndex = i;
      int toVectorNo = 0;
      int toVectorIndex = 0;
      bool avoidCopyIfPossible = true;
      bool slotIsConnected = getSlotInformation(fromVectorNo, fromVectorIndex, toVectorNo, toVectorIndex, avoidCopyIfPossible);

      if (!slotIsConnected)
        continue;

      int fromIndex = fromVectorNo*nFieldVariablesTerm2Vector1_ + fromVectorIndex;
      connectorTerm2To1_[fromIndex].avoidCopyIfPossible = false;
    }
  }

  LOG(DEBUG) << "now initialize lookup-table.";

  // fill look-up table slotInformation_[transferDirectionTerm1To2][fromVectorNo][fromVectorIndex]
  for (int transferDirectionTerm1To2 = 0; transferDirectionTerm1To2 < 2; transferDirectionTerm1To2++)
  {
    for (int fromVectorNo = 0; fromVectorNo < 2; fromVectorNo++)
    {
      // iterate to the next fromVectorIndex until there were 3 unconnected slots, then it is assumed that there are no more
      int nUnsuccessful = 0;
      for (int fromVectorIndex = 0; ; fromVectorIndex++)
      {
        transferDirectionTerm1To2_ = transferDirectionTerm1To2;
        Result result;
        result.successful = getSlotInformation(fromVectorNo, fromVectorIndex, result.toVectorNo, result.toVectorIndex, result.avoidCopyIfPossible, true);

        slotInformation_[transferDirectionTerm1To2][fromVectorNo].push_back(result);

        if (!result.successful)
          nUnsuccessful++;

        if (nUnsuccessful == 3)
          break;
      }
    }
  }

  slotInformationInitialized_ = true;
}
