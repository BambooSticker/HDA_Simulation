## Overview
The HDA Hydrodealkylation process is an important way to produce benzene. 
The toluene presented to excess hydrogen is dealkylated to form benzene 
with the reaction: 
$$C_6H_5CH_3 + H_2 \rightarrow C_6H_6 + CH_4$$
<p>
  Learning on the <a href="https://github.com/IDAES/idaes-pse" target="_blank">IDAES project</a>, 
  in this repository, a simplified HDA process is constructed, simulated and optimized using Pyomo. 
  More tools (pandas, matplotlib, etc.) are leveraged to exploit and visualize the simulation outcome. 
</p>

## Introduction
A brief introduction of this project is as follows:
-  **hda_flowsheet.py** and **HDA_Flowsheet.ipynb** contain the complete simulation and optimization process. A ConcreteModel is
constructed as a blank flowsheet, then every unit model involved is created and linked with Arc. After
selecting tear stream and tearing the cycle, Sequential Decomposition was leveraged to initialize the model.
When initial simulation is successful, the operation conditions of several units are relaxed as decision variables,
and operation cost based on utility usage can be defined as the objective of the large-scale optimization problem. 
- **State_Block.py** is used to create some thermophysical properties for units, including Antoine equations,
bubble and dew approximation, smooth VLE calculation, etc. Due to limited coding skills and time (XD), the state
of Mixer unit is created using another **State_Block_for_Mixer.py** for its more-than-one inlets.
- **Unit_Model.py** is a library of operation units involved in the HDA process, including flash, heater, mixer,
spliter and reactor. A unit model consists of a state block that contains thermophysical properties, and a set of
balance constraints depending on the spesific type of unit.
- In **Method.py**, several customized and special strategies are defined as functions to contribute to
successful initialization and visualization. 

## About
<p>
  This repository is for the undergraduate final thesis at <strong>East China University of Science and Technology</strong>. 
</p>
</div>


