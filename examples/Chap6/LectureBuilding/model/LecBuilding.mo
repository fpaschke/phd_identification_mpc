within ;
package LecBuilding "Several Utilities"
  block NPersToHeatFlow
    "Number of Persons in HeatFlow inputs to thermal Zone"
    extends Modelica.Blocks.Icons.Block;
    parameter Real[4,1] A = [1; 1; 1; 1] "Area of Hs1 to Hs4 m^2";
    parameter Real[3,1] Q = [70; 15; 15] "Radiant covective and latent Heat Gain per Person in W";

    Modelica.Blocks.Interfaces.RealInput u[4,1]
      "4x1 Real input signal" annotation (Placement(transformation(
            extent={{-138,-20},{-98,20}})));
    Modelica.Blocks.Interfaces.RealOutput y1[3,1]
      "Connector of Real output signals 1" annotation (Placement(transformation(
            extent={{100,80},{120,100}})));
    Modelica.Blocks.Interfaces.RealOutput y2[3,1]
      "Connector of Real output signals 2" annotation (Placement(transformation(
            extent={{100,20},{120,40}})));
    Modelica.Blocks.Interfaces.RealOutput y3[3,1]
      "Connector of Real output signals 3" annotation (Placement(transformation(
            extent={{100,-40},{120,-20}})));
    Modelica.Blocks.Interfaces.RealOutput y4[3,1]
      "Connector of Real output signals 4" annotation (Placement(transformation(
            extent={{100,-100},{120,-80}})));
  equation
    y1 = Q*u[1,1]/A[1,1];
    y2 = Q*u[2,1]/A[2,1];
    y3 = Q*u[3,1]/A[3,1];
    y4 = Q*u[4,1]/A[4,1];
  end NPersToHeatFlow;

  model TestLecBuilding "Test of LectureBuilding Model"
    LectureBuilding LecBuil
      annotation (Placement(transformation(extent={{74,-12},{96,8}})));
    Modelica.Blocks.Routing.Multiplex4 multiplex4_1
      annotation (Placement(transformation(extent={{-12,20},{8,40}})));
    Modelica.Blocks.Math.MatrixGain gai(K=0.5*[600; 600; 300; 300])
      "Number of Persons"
      annotation (Placement(transformation(extent={{40,-28},{60,-8}})));
    Modelica.Blocks.Sources.Pulse nPer(
      width=10/24*100,
      amplitude=1,
      startTime(displayUnit="h") = 28800,
      period(displayUnit="h") = 86400) "Number of persons"
      annotation (Placement(transformation(extent={{-6,-28},{6,-16}})));
    Modelica.Blocks.Sources.RealExpression TsupAHU(y=19)
      annotation (Placement(transformation(extent={{-12,-6},{8,14}})));
    Modelica.Blocks.Sources.Pulse m1(
      width=13/24*100,
      amplitude=6,
      startTime(displayUnit="h") = 21600,
      period(displayUnit="h") = 86400) "Number of persons"
      annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
    Modelica.Blocks.Sources.Pulse m2(
      width=13/24*100,
      amplitude=6,
      startTime(displayUnit="h") = 21600,
      period(displayUnit="h") = 86400) "Number of persons"
      annotation (Placement(transformation(extent={{-80,10},{-60,30}})));
    Modelica.Blocks.Sources.Pulse m3(
      width=13/24*100,
      amplitude=6,
      startTime(displayUnit="h") = 21600,
      period(displayUnit="h") = 86400) "Number of persons"
      annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));
    Modelica.Blocks.Sources.Pulse m4(
      width=13/24*100,
      amplitude=6,
      startTime(displayUnit="h") = 21600,
      period(displayUnit="h") = 86400) "Number of persons"
      annotation (Placement(transformation(extent={{-80,-54},{-60,-34}})));
    Modelica.Blocks.Sources.RealExpression Qhall(y=0)
      annotation (Placement(transformation(extent={{16,-78},{36,-58}})));
    Modelica.Blocks.Sources.Pulse wd(
      width=5/7*100,
      amplitude=1,
      startTime(displayUnit="d") = 172800,
      period(displayUnit="d") = 604800) "Weekday/Weekend"
      annotation (Placement(transformation(extent={{-80,74},{-60,94}})));
    Modelica.Blocks.Math.Product product
      annotation (Placement(transformation(extent={{-46,42},{-36,52}})));
    Modelica.Blocks.Math.Product product1
      annotation (Placement(transformation(extent={{-46,12},{-36,22}})));
    Modelica.Blocks.Math.Product product2
      annotation (Placement(transformation(extent={{-46,-18},{-36,-8}})));
    Modelica.Blocks.Math.Product product3
      annotation (Placement(transformation(extent={{-46,-52},{-36,-42}})));
    Modelica.Blocks.Math.Product product4
      annotation (Placement(transformation(extent={{20,-24},{32,-12}})));
  equation
    connect(multiplex4_1.y,LecBuil.M_Sup)  annotation (Line(points={{9,30},{34,
            30},{34,5.5},{73.0435,5.5}}, color={0,0,127}));
    connect(gai.y, LecBuil.N_Pers) annotation (Line(points={{61,-18},{64,-18},{
            64,-2.5},{73.0435,-2.5}}, color={0,0,127}));
    connect(LecBuil.T_Sup, TsupAHU.y)
      annotation (Line(points={{73.0435,3.7},{9,3.7},{9,4}}, color={0,0,127}));
    connect(Qhall.y, LecBuil.Q_Hall) annotation (Line(points={{37,-68},{68,-68},
            {68,-10},{73.0435,-10}}, color={0,0,127}));
    connect(m1.y, product.u1)
      annotation (Line(points={{-59,50},{-47,50}}, color={0,0,127}));
    connect(m2.y, product1.u1)
      annotation (Line(points={{-59,20},{-47,20}}, color={0,0,127}));
    connect(m3.y, product2.u1)
      annotation (Line(points={{-59,-10},{-47,-10}}, color={0,0,127}));
    connect(m4.y, product3.u1) annotation (Line(points={{-59,-44},{-52,-44},{
            -52,-44},{-47,-44}}, color={0,0,127}));
    connect(product3.y, multiplex4_1.u4[1]) annotation (Line(points={{-35.5,-47},
            {-22,-47},{-22,21},{-14,21}}, color={0,0,127}));
    connect(product2.y, multiplex4_1.u3[1]) annotation (Line(points={{-35.5,-13},
            {-26,-13},{-26,27},{-14,27}}, color={0,0,127}));
    connect(product1.y, multiplex4_1.u2[1]) annotation (Line(points={{-35.5,17},
            {-28,17},{-28,33},{-14,33}}, color={0,0,127}));
    connect(product.y, multiplex4_1.u1[1]) annotation (Line(points={{-35.5,47},
            {-26,47},{-26,39},{-14,39}}, color={0,0,127}));
    connect(wd.y, product.u2) annotation (Line(points={{-59,84},{-52,84},{-52,
            44},{-47,44}}, color={0,0,127}));
    connect(product1.u2, product.u2) annotation (Line(points={{-47,14},{-52,14},
            {-52,44},{-47,44}}, color={0,0,127}));
    connect(product2.u2, product.u2) annotation (Line(points={{-47,-16},{-52,
            -16},{-52,44},{-47,44}}, color={0,0,127}));
    connect(product3.u2, product.u2) annotation (Line(points={{-47,-50},{-52,
            -50},{-52,44},{-47,44}}, color={0,0,127}));
    connect(nPer.y, product4.u2) annotation (Line(points={{6.6,-22},{8,-22},{8,
            -21.6},{18.8,-21.6}}, color={0,0,127}));
    connect(product4.y, gai.u[1])
      annotation (Line(points={{32.6,-18},{38,-18}}, color={0,0,127}));
    connect(product4.u1, product.u2) annotation (Line(points={{18.8,-14.4},{16,
            -14.4},{16,84},{-52,84},{-52,44},{-47,44}}, color={0,0,127}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StopTime=31536000,
        Interval=600,
        __Dymola_fixedstepsize=600,
        __Dymola_Algorithm="Dassl"));
  end TestLecBuilding;

  model MultiLayer
    "Copied from Buildings Library: Changed intialization in SingleLayer"
    extends Buildings.HeatTransfer.Conduction.BaseClasses.PartialConductor(
     final R=sum(lay[i].R for i in 1:nLay));
    Modelica.SIunits.Temperature T[sum(layers.nSta)](
      each nominal = 300) "Temperature at the states";
    Modelica.SIunits.HeatFlowRate Q_flow[sum(layers.nSta)+nLay]
      "Heat flow rate from state i to i+1";
    extends Buildings.HeatTransfer.Conduction.BaseClasses.PartialConstruction;

    parameter Boolean stateAtSurface_a=true
      "=true, a state will be at the surface a"
      annotation (Dialog(tab="Dynamics"),
                  Evaluate=true);
    parameter Boolean stateAtSurface_b=true
      "=true, a state will be at the surface b"
      annotation (Dialog(tab="Dynamics"),
                  Evaluate=true);

  protected
    LecBuilding.SingleLayer[nLay] lay(
     final nSta2={layers.nSta[i] for i in 1:nLay},
     each final A=A,
     final stateAtSurface_a = {if i == 1 then stateAtSurface_a else false for i in 1:nLay},
     final stateAtSurface_b = {if i == nLay then stateAtSurface_b else false for i in 1:nLay},
     material = {layers.material[i] for i in 1:size(layers.material, 1)},
     T_a_start = { T_b_start+(T_a_start-T_b_start) * 1/R *
      sum(layers.material[k].R for k in i:size(layers.material, 1)) for i in 1:size(layers.material, 1)},
     T_b_start = { T_a_start+(T_b_start-T_a_start) * 1/R *
      sum(layers.material[k].R for k in 1:i) for i in 1:size(layers.material, 1)},
     each steadyStateInitial = steadyStateInitial) "Material layer"
      annotation (Placement(transformation(extent={{-20,-10},{0,10}})));

  equation
    // This section assigns the temperatures and heat flow rates of the layer models to
    // an array that makes plotting the results easier.
    for i in 1:nLay loop
      for j in 1:layers.nSta[i] loop
        T[sum(layers.nSta[k] for k in 1:(i-1)) +j] = lay[i].T[j];
      end for;
      for j in 1:layers.nSta[i]+1 loop
        Q_flow[sum(layers.nSta[k] for k in 1:i-1)+(i-1)+j] = lay[i].Q_flow[j];
      end for;
    end for;
    connect(port_a, lay[1].port_a) annotation (Line(
        points={{-100,5.55112e-16},{-60,5.55112e-16},{-60,6.10623e-16},{-20,
            6.10623e-16}},
        color={191,0,0},
        smooth=Smooth.None));
    for i in 1:nLay-1 loop
    connect(lay[i].port_b, lay[i+1].port_a) annotation (Line(
        points={{5.55112e-16,6.10623e-16},{20,6.10623e-16},{20,-20},{-40,-20},{
              -40,6.10623e-16},{-20,6.10623e-16}},
        color={191,0,0},
        smooth=Smooth.None));
    end for;
    connect(lay[nLay].port_b, port_b) annotation (Line(
        points={{5.55112e-16,6.10623e-16},{49,6.10623e-16},{49,5.55112e-16},{100,
            5.55112e-16}},
        color={191,0,0},
        smooth=Smooth.None));

    annotation ( Icon(coordinateSystem(
            preserveAspectRatio=true, extent={{-100,-100},{100,100}}), graphics={
     Rectangle(
      extent={{0,80},{80,-80}},       fillColor={175,175,175},
     fillPattern=FillPattern.Solid,    lineColor={175,175,175}),
     Rectangle(
      extent={{-80,80},{0,-80}},      fillColor={215,215,215},
     fillPattern=FillPattern.Solid,    lineColor={175,175,175}),
     Line(points={{-92,0},{90,0}},      color = {0, 0, 0}, thickness = 0.5,
     smooth = Smooth.None),
     Line(points={{-18,-40},{-32,-40}},     color = {0, 0, 0}, thickness = 0.5,
     smooth = Smooth.None),
     Line(points={{-12,-32},{-38,-32}},     color = {0, 0, 0}, thickness = 0.5,
     smooth = Smooth.None),            Line(points={{-25,0},{-25,-32}},
     color = {0, 0, 0}, thickness = 0.5, smooth = Smooth.None),
     Line(points={{32,-40},{18,-40}},       color = {0, 0, 0}, thickness = 0.5,
     smooth = Smooth.None),
     Line(points={{38,-32},{12,-32}},       color = {0, 0, 0}, thickness = 0.5,
     smooth = Smooth.None),            Line(points={{25,0},{25,-32}},
     color = {0, 0, 0}, thickness = 0.5, smooth = Smooth.None),
                                       Rectangle(extent={{-60,6},{-40,-6}},
     lineColor = {0, 0, 0}, lineThickness =  0.5, fillColor = {255, 255, 255},
     fillPattern = FillPattern.Solid), Rectangle(extent={{-10,6},{10,-6}},
     lineColor = {0, 0, 0}, lineThickness =  0.5, fillColor = {255, 255, 255},
     fillPattern = FillPattern.Solid), Rectangle(extent={{40,6},{60,-6}},
     lineColor = {0, 0, 0}, lineThickness =  0.5, fillColor = {255, 255, 255},
     fillPattern = FillPattern.Solid),
     Line(points={{86,-40},{72,-40}},       color = {0, 0, 0}, thickness = 0.5,
     smooth = Smooth.None,
     visible=stateAtSurface_b),
     Line(points={{92,-32},{66,-32}},       color = {0, 0, 0}, thickness = 0.5,
     smooth = Smooth.None,
     visible=stateAtSurface_b),            Line(points={{79,0},{79,-32}},
     color = {0, 0, 0}, thickness = 0.5, smooth = Smooth.None,
     visible=stateAtSurface_b),
     Line(points={{-79,0},{-79,-32}},
     color = {0, 0, 0}, thickness = 0.5, smooth = Smooth.None,
     visible=stateAtSurface_a),
     Line(points={{-66,-32},{-92,-32}},     color = {0, 0, 0}, thickness = 0.5,
     smooth = Smooth.None,
     visible=stateAtSurface_a),
     Line(points={{-72,-40},{-86,-40}},     color = {0, 0, 0}, thickness = 0.5,
     smooth = Smooth.None,
     visible=stateAtSurface_a)}),
      defaultComponentName="heaCon",
      Documentation(info="<html>
<p>
This is a model of a heat conductor with multiple material layers and energy storage.
The construction has at least one material layer, and each layer has
at least one temperature node. The layers are modeled using an instance of
<a href=\"Buildings.HeatTransfer.Conduction.SingleLayer\">
Buildings.HeatTransfer.Conduction.SingleLayer</a>.
See this model for an explanation of the equations that are applied to
each material layer.
</p>
<h4>Important parameters</h4>
<p>
The construction material is defined by a record of the package
<a href=\"modelica://Buildings.HeatTransfer.Data.OpaqueConstructions\">
Buildings.HeatTransfer.Data.OpaqueConstructions</a>.
This record allows specifying materials that store energy, and material
that are a thermal conductor only with no heat storage.
To assign the material properties to this model, do the following:
</p>
<ol>
<li>
Create an instance of a record of
<a href=\"modelica://Buildings.HeatTransfer.Data.OpaqueConstructions\">
Buildings.HeatTransfer.Data.OpaqueConstructions</a>, for example
by dragging the record into the schematic model editor.
</li>
<li>
Make sure the instance has the attribute <code>parameter</code>, which may not be
assigned automatically when you drop the model in a graphical editor. For
example, an instanciation may look like
<pre>
 parameter Data.OpaqueConstructions.Insulation100Concrete200 layers
   \"Material layers of construction\"
   annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
</pre>
</li>
<li>
Assign the instance of the material to the instance of the heat transfer
model as shown in
<a href=\"modelica://Buildings.HeatTransfer.Examples.ConductorMultiLayer\">
Buildings.HeatTransfer.Examples.ConductorMultiLayer</a>.
</li>
</ol>
<p>
The parameters <code>stateAtSurface_a</code> and
<code>stateAtSurface_b</code>
determine whether there is a state variable at these surfaces,
as described above.
Note that if <code>stateAtSurface_a = true</code>,
then there is temperature state on the surface a with prescribed
value, as determined by the differential equation of the heat conduction.
Hence, in this situation, it is not possible to
connect a temperature boundary condition such as
<a href=\"modelica://Buildings.HeatTransfer.Sources.FixedTemperature\">
Buildings.HeatTransfer.Sources.FixedTemperature</a> as this would
yield to specifying the same temperature twice.
To avoid this, either set <code>stateAtSurface_a = false</code>,
or place a thermal resistance
between the boundary condition and the surface of this model.
The same applies for surface b.
See the examples in
<a href=\"modelica://Buildings.HeatTransfer.Examples\">
Buildings.HeatTransfer.Examples</a>.
</p>
</html>",   revisions="<html>
<ul>
<li>
October 16, 2017, by Michael Wetter:<br/>
Corrected wrong result variable <code>R</code> and <code>UA</code>.
These variables are only used for reporting.
All other calculations were not affected by this error.
</li>
<li>
January 05, 2017, by Thierry S. Nouidui:<br/>
Removed parameter <code>nSta2</code>.
</li>
<li>
November 17, 2016, by Thierry S. Nouidui:<br/>
Added parameter <code>nSta2</code> to avoid translation error
in Dymola 2107. This is a work-around for a bug in Dymola
which will be addressed in future releases.
</li>
<li>
October 29, 2016, by Michael Wetter:<br/>
Added option to place a state at the surface.<br/>
This is for
<a href=\"https://github.com/lbl-srg/modelica-buildings/issues/565\">issue 565</a>.
</li>
<li>
September 24, 2015 by Michael Wetter:<br/>
Set the start value of <code>T</code>.
This is for
<a href=\"https://github.com/lbl-srg/modelica-buildings/issues/426\">issue 426</a>.
</li>
<li>
March 18, 2015, by Michael Wetter:<br/>
Replaced <code>nLay</code> in the <code>sum()</code> of the parameter assignment
with <code>size(layers.material, 1)</code> to avoid incorrect results in OpenModelica.
See <a href=\"https://github.com/lbl-srg/modelica-buildings/commit/4578a3d3b80e760cc83d705963f3b17e41c1e7da#diff-9628c0eecd08caed8b30f1f993de7501L12\">github note</a>.
</li>
<li>
March 13, 2015, by Michael Wetter:<br/>
Changed assignment of <code>nLay</code> to avoid a translation error
in OpenModelica.
</li>
<li>
October 15, 2014, by Michael Wetter:<br/>
Changed assignment of <code>R</code> to be in the <code>extends</code> statement
to avoid a division by zero in OpenModelica.
</li>
<li>
September 9, 2014, by Michael Wetter:<br/>
Reverted change from March 1 2013 as this causes an error during model check
in Dymola 2015 FD01 beta1.
</li>
<li>
August 12, 2014, by Michael Wetter:<br/>
Reformulated the protected elements and the model instantiation to avoid
a warning in the OpenModelica parser.
</li>
<li>
March 1, 2013, by Michael Wetter:<br/>
Removed <code>initial equation</code> section and assigned the protected parameters
<code>_T_a_start</code> and <code>_T_b_start</code> directly to avoid a warning during
translation.
</li>
<li>
March 6 2010, by Michael Wetter:<br/>
First implementation.
</li>
</ul>
</html>"));
  end MultiLayer;

  model SingleLayer "Model for single layer heat conductance"
    extends Buildings.HeatTransfer.Conduction.BaseClasses.PartialConductor(
     final R=if (material.R < Modelica.Constants.eps) then material.x/material.k/A else material.R/A);
     // if material.R == 0, then the material specifies material.k, and this model specifies x
     // For resistances, material.k need not be specified, and hence we use material.R

    // The value T[:].start is used by the solver when finding initial states
    // that satisfy dT/dt=0, which requires solving a system of nonlinear equations
    // if the convection coefficient is a function of temperature.
    Modelica.SIunits.Temperature T[nSta](start=
     if stateAtSurface_a then
       cat(1,
         {T_a_start},
         {(T_a_start + (T_b_start - T_a_start)*UA*sum(RNod[k] for k in 1:i-1)) for i in 2:nSta})
     else
      {(T_a_start + (T_b_start - T_a_start)*UA*sum(RNod[k] for k in 1:i)) for i in 1:nSta},
     each nominal=300)
      "Temperature at the states";

    Modelica.SIunits.HeatFlowRate Q_flow[nSta+1](each start=0)
      "Heat flow rates to each state";
    Modelica.SIunits.SpecificInternalEnergy u[nSta](
      each start=2.7E5,
      each nominal=2.7E5)
      "Definition of specific internal energy";

    parameter Boolean stateAtSurface_a=true
      "=true, a state will be at the surface a"
      annotation (Dialog(tab="Dynamics"),
                  Evaluate=true);
    parameter Boolean stateAtSurface_b=true
      "=true, a state will be at the surface b"
      annotation (Dialog(tab="Dynamics"),
                  Evaluate=true);

    replaceable parameter Buildings.HeatTransfer.Data.BaseClasses.Material material
      "Material from Data.Solids, Data.SolidsPCM or Data.Resistances" annotation (
       choicesAllMatching=true, Placement(transformation(extent={{60,60},{80,80}})));

    parameter Boolean steadyStateInitial=false
      "=true initializes dT(0)/dt=0, false initializes T(0) at fixed temperature using T_a_start and T_b_start"
          annotation (Dialog(group="Initialization"), Evaluate=true);
    parameter Modelica.SIunits.Temperature T_a_start=293.15
      "Initial temperature at port_a, used if steadyStateInitial = false"
      annotation (Dialog(group="Initialization", enable=not steadyStateInitial));
    parameter Modelica.SIunits.Temperature T_b_start=293.15
      "Initial temperature at port_b, used if steadyStateInitial = false"
      annotation (Dialog(group="Initialization", enable=not steadyStateInitial));
    parameter Integer nSta2=material.nSta
    "Number of states in a material (do not overwrite, used to work around Dymola 2017 bug)"
       annotation (Evaluate=true, HideResult=true, Dialog(enable=false, tab="Advanced"));
  protected
    final parameter Integer nSta=
      max(nSta2,
          if stateAtSurface_a or stateAtSurface_b then 2 else 1)
      "Number of state variables";
    final parameter Integer nR=nSta+1 "Number of thermal resistances";
    parameter Modelica.SIunits.ThermalResistance RNod[nR]=
      if (stateAtSurface_a and stateAtSurface_b) then
        if (nSta==2) then
          {(if i==1 or i==nR then 0 else R/(nSta-1)) for i in 1:nR}
        else
          {(if i==1 or i==nR then 0 elseif i==2 or i==nR-1 then R/(2*(nSta-2)) else R/(nSta-2)) for i in 1:nR}
        elseif (stateAtSurface_a and (not stateAtSurface_b)) then
          {(if i==1 then 0 elseif i==2 or i==nR then R/(2*(nSta-1)) else R/(nSta-1)) for i in 1:nR}
      elseif (stateAtSurface_b and (not stateAtSurface_a)) then
         {(if i==nR then 0 elseif i==1 or i==nR-1 then R/(2*(nSta-1)) else R/(nSta-1)) for i in 1:nR}
      else
        {R/(if i==1 or i==nR then (2*nSta) else nSta) for i in 1:nR}
      "Thermal resistance";

    parameter Modelica.SIunits.Mass m[nSta]=
     (A*material.x*material.d) *
     (if (stateAtSurface_a and stateAtSurface_b) then
       if (nSta==2) then
         {1/(2*(nSta-1)) for i in 1:nSta}
       elseif (nSta==3) then
         {1/(if i==1 or i==nSta then (2*(nSta-1)) else (nSta-1)) for i in 1:nSta}
       else
         {1/(if i==1 or i==nSta or i==2 or i==nSta-1 then (2*(nSta-2)) else (nSta-2)) for i in 1:nSta}
       elseif (stateAtSurface_a and (not stateAtSurface_b)) then
         {1/(if i==1 or i==2 then (2*(nSta-1)) else (nSta-1)) for i in 1:nSta}
       elseif (stateAtSurface_b and (not stateAtSurface_a)) then
         {1/(if i==nSta or i==nSta-1 then (2*(nSta-1)) else (nSta-1)) for i in 1:nSta}
       else
         {1/(nSta) for i in 1:nSta})
      "Mass associated with the temperature state";

    final parameter Real mInv[nSta]=
      if material.steadyState then zeros(nSta) else {1/m[i] for i in 1:nSta}
      "Inverse of the mass associated with the temperature state";

    final parameter Modelica.SIunits.HeatCapacity C[nSta] = m*material.c
      "Heat capacity associated with the temperature state";
    final parameter Real CInv[nSta]=
      if material.steadyState then zeros(nSta) else {1/C[i] for i in 1:nSta}
      "Inverse of heat capacity associated with the temperature state";

    parameter Modelica.SIunits.SpecificInternalEnergy ud[Buildings.HeatTransfer.Conduction.nSupPCM](
      each fixed=false)
      "Support points for derivatives (used for PCM)";
    parameter Modelica.SIunits.Temperature Td[Buildings.HeatTransfer.Conduction.nSupPCM](
      each fixed=false)
      "Support points for derivatives (used for PCM)";
    parameter Real dT_du[Buildings.HeatTransfer.Conduction.nSupPCM](
      each fixed=false,
      each unit="kg.K2/J")
      "Derivatives dT/du at the support points (used for PCM)";

  initial equation
    // The initialization is only done for materials that store energy.
      if not material.steadyState then
        if steadyStateInitial then
          if material.phasechange then
            der(u) = zeros(nSta);
          else
            der(T) = zeros(nSta);
          end if;
        else
          if stateAtSurface_a then
            // T[1] = T_a_start;
            // for i in 2:nSta loop
            for i in 2:(nSta-1) loop
              T[i] =T_a_start + (T_b_start - T_a_start)*UA*sum(RNod[k] for k in 1:i-1);
            end for;
          else // stateAtSurface_a == false
            for i in 1:nSta loop
              T[i] = T_a_start + (T_b_start - T_a_start)*UA*sum(RNod[k] for k in 1:i);
            end for;
          end if;
        end if;
      end if;

     if material.phasechange then
       (ud, Td, dT_du) = Buildings.HeatTransfer.Conduction.BaseClasses.der_temperature_u(
         c =  material.c,
         TSol=material.TSol,
         TLiq=material.TLiq,
         LHea=material.LHea,
         ensureMonotonicity=material.ensureMonotonicity);
     else
       ud    = zeros(Buildings.HeatTransfer.Conduction.nSupPCM);
       Td    = zeros(Buildings.HeatTransfer.Conduction.nSupPCM);
       dT_du = zeros(Buildings.HeatTransfer.Conduction.nSupPCM);
     end if;
  equation
      port_a.Q_flow = +Q_flow[1];
      port_b.Q_flow = -Q_flow[end];

      port_a.T-T[1]    = if stateAtSurface_a then 0 else Q_flow[1]*RNod[1];
      T[nSta]-port_b.T = if stateAtSurface_b then 0 else Q_flow[end]*RNod[end];

      for i in 1:nSta-1 loop
         // Q_flow[i+1] is heat flowing from (i) to (i+1)
         // because T[1] has Q_flow[1] and Q_flow[2] acting on it.
         T[i]-T[i+1] = Q_flow[i+1]*RNod[i+1];
      end for;

      // Steady-state heat balance
      if material.steadyState then
        for i in 2:nSta+1 loop
          Q_flow[i] = port_a.Q_flow;
        end for;

        for i in 1:nSta loop
          if material.phasechange then
            // Phase change material
            T[i]=Buildings.HeatTransfer.Conduction.BaseClasses.temperature_u(
                      ud=ud,
                      Td=Td,
                      dT_du=dT_du,
                      u=u[i]);
          else
            // Regular material
            u[i]=0; // u is not required in this case
          end if;
        end for;
      else
        // Transient heat conduction
        if material.phasechange then
          // Phase change material
          for i in 1:nSta loop
            der(u[i]) = (Q_flow[i]-Q_flow[i+1])*mInv[i];
            // Recalculation of temperature based on specific internal energy
            T[i]=Buildings.HeatTransfer.Conduction.BaseClasses.temperature_u(
                      ud=ud,
                      Td=Td,
                      dT_du=dT_du,
                      u=u[i]);
          end for;
        else
          // Regular material
          for i in 1:nSta loop
            der(T[i]) = (Q_flow[i]-Q_flow[i+1])*CInv[i];
          end for;
          for i in 1:nSta loop
            u[i]=0; // u is not required in this case
          end for;
        end if;
      end if;

    annotation ( Icon(coordinateSystem(
            preserveAspectRatio=false,extent={{-100,-100},{100,100}}), graphics={
          Text(
            extent={{-100,-80},{6,-98}},
            lineColor={0,0,255},
            textString="%material.x"),
          Text(
            extent={{8,-74},{86,-104}},
            lineColor={0,0,255},
            textString="%nSta"),
     Rectangle(
      extent={{-60,80},{60,-80}},     fillColor={215,215,215},
     fillPattern=FillPattern.Solid,    lineColor={175,175,175}),
     Line(points={{-92,0},{90,0}},      color = {0, 0, 0}, thickness = 0.5,
     smooth = Smooth.None),
     Line(points={{8,-40},{-6,-40}},        color = {0, 0, 0}, thickness = 0.5,
     smooth = Smooth.None),
     Line(points={{14,-32},{-12,-32}},      color = {0, 0, 0}, thickness = 0.5,
     smooth = Smooth.None),            Line(
            points={{0,0},{0,-32}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None),       Rectangle(extent={{-40,6},{-20,-6}},
     lineColor = {0, 0, 0}, lineThickness =  0.5, fillColor = {255, 255, 255},
     fillPattern = FillPattern.Solid), Rectangle(extent={{20,6},{40,-6}},
     lineColor = {0, 0, 0}, lineThickness =  0.5, fillColor = {255, 255, 255},
     fillPattern = FillPattern.Solid),
     Line(points={{66,-40},{52,-40}},       color = {0, 0, 0}, thickness = 0.5,
     smooth = Smooth.None,
     visible=stateAtSurface_b),
     Line(points={{72,-32},{46,-32}},       color = {0, 0, 0}, thickness = 0.5,
     smooth = Smooth.None,
     visible=stateAtSurface_b),            Line(points={{59,0},{59,-32}},
     color = {0, 0, 0}, thickness = 0.5, smooth = Smooth.None,
     visible=stateAtSurface_b),
     Line(points={{-59,0},{-59,-32}},
     color = {0, 0, 0}, thickness = 0.5, smooth = Smooth.None,
     visible=stateAtSurface_a),
     Line(points={{-46,-32},{-72,-32}},     color = {0, 0, 0}, thickness = 0.5,
     smooth = Smooth.None,
     visible=stateAtSurface_a),
     Line(points={{-52,-40},{-66,-40}},     color = {0, 0, 0}, thickness = 0.5,
     smooth = Smooth.None,
     visible=stateAtSurface_a)}),
  defaultComponentName="lay",
      Documentation(info="<html>
<p>
This is a model of a heat conductor for a single layer of homogeneous material
that computes transient or steady-state heat conduction.
</p>
<h4>Main equations</h4>
<h5>Transient heat conduction in materials without phase change</h5>
<p>
If the material is a record that extends
<a href=\"modelica://Buildings.HeatTransfer.Data.Solids\">
Buildings.HeatTransfer.Data.Solids</a> and its
specific heat capacity (as defined by the record <code>material.c</code>)
is non-zero, then this model computes <i>transient</i> heat conduction, i.e., it
computes a numerical approximation to the solution of the heat equation
</p>
<p align=\"center\" style=\"font-style:italic;\">
   &rho; c (&part; T(s,t) &frasl; &part;t) =
    k (&part;&sup2; T(s,t) &frasl; &part;s&sup2;),
</p>
<p>
where
<i>&rho;</i>
is the mass density,
<i>c</i>
is the specific heat capacity per unit mass,
<i>T</i>
is the temperature at location <i>s</i> and time <i>t</i> and
<i>k</i> is the heat conductivity.
At the locations <i>s=0</i> and <i>s=x</i>, where <i>x</i> is the
material thickness, the temperature and heat flow rate is equal to the
temperature and heat flow rate of the heat ports.
</p>
<h5>Transient heat conduction in phase change materials</h5>
<p>
If the material is declared using a record of type
<a href=\"modelica://Buildings.HeatTransfer.Data.SolidsPCM\">
Buildings.HeatTransfer.Data.SolidsPCM</a>, the heat transfer
in a phase change material is computed.
The record <a href=\"modelica://Buildings.HeatTransfer.Data.SolidsPCM\">
Buildings.HeatTransfer.Data.SolidsPCM</a>
declares the solidus temperature <code>TSol</code>,
the liquidus temperature <code>TLiq</code> and the latent heat of
phase transformation <code>LHea</code>.
For heat transfer with phase change, the specific internal energy <i>u</i>
is the dependent variable, rather than the temperature.
Therefore, the governing equation is
</p>
<p align=\"center\" style=\"font-style:italic;\">
   &rho; (&part; u(s,t) &frasl; &part;t) =
    k (&part;&sup2; T(s,t) &frasl; &part;s&sup2;).
</p>
<p>
The constitutive
relation between specific internal energy <i>u</i> and temperature <i>T</i> is defined in
<a href=\"modelica://Buildings.HeatTransfer.Conduction.BaseClasses.temperature_u\">
Buildings.HeatTransfer.Conduction.BaseClasses.temperature_u</a> by using
cubic hermite spline interpolation with linear extrapolation.
</p>
<h5>Steady-state heat conduction</h5>
<p>
If <code>material.c=0</code>, or if the material extends
<a href=\"modelica://Buildings.HeatTransfer.Data.Resistances\">
Buildings.HeatTransfer.Data.Resistances</a>,
then steady-state heat conduction is computed. In this situation, the heat
flow between its heat ports is
</p>
<p align=\"center\" style=\"font-style:italic;\">
   Q = A &nbsp; k &frasl; x &nbsp; (T<sub>a</sub>-T<sub>b</sub>),
</p>
<p>
where
<i>A</i> is the cross sectional area,
<i>x</i> is the layer thickness,
<i>T<sub>a</sub></i> is the temperature at port a and
<i>T<sub>b</sub></i> is the temperature at port b.
</p>
<h5>Spatial discretization</h5>
<p>
To spatially discretize the heat equation, the construction is
divided into compartments (control volumes) with <code>material.nSta &ge; 1</code> state variables.
Each control volume has the same material properties.
The state variables are connected to each other through thermal resistances.
If <code>stateAtSurface_a = true</code>, a state is placed
at the surface a, and similarly, if
<code>stateAtSurface_b = true</code>, a state is placed
at the surface b.
Otherwise, these states are placed inside the material, away
from the surface.
Thus, to obtain
the surface temperature, use <code>port_a.T</code> (or <code>port_b.T</code>)
and not the variable <code>T[1]</code>.
</p>

As an example, we assume a material with a length of <code>x</code>
and a discretization with four state variables.
<ul>
<li>
If <code>stateAtSurface_a = false</code> and <code>stateAtSurface_b = false</code>,
then each of the four state variables is placed in the middle of a control volume with length <code>l=x/material.nSta</code>.
<p align=\"left\"><img alt=\"image\" src=\"modelica://Buildings/Resources/Images/HeatTransfer/Conduction/noStateAtSurface.svg\"/>
</li>
<li>
If <code>stateAtSurface_a = true</code> or <code>stateAtSurface_b = true</code>,
then one state is placed on the surface of the material. Each of the remaining three states
is placed in the middle of a control volume with length <code>l=x/(material.nSta-1)</code>.
<p align=\"left\"><img alt=\"image\" src=\"modelica://Buildings/Resources/Images/HeatTransfer/Conduction/oneStateAtSurface.svg\"/>
</li>
<li>
If <code>stateAtSurface_a = true</code> and <code>stateAtSurface_b = true</code>,
then two states are placed on the surfaces of the material. Each of the remaining two states is placed
in the middle of a control volume with length <code>l=x/(material.nSta-2)</code>.
<p align=\"left\"><img alt=\"image\" src=\"modelica://Buildings/Resources/Images/HeatTransfer/Conduction/twoStatesAtSurface.svg\"/>
</li>
</ul>

<p>
To build multi-layer constructions,
use
<a href=\"Buildings.HeatTransfer.Conduction.MultiLayer\">
Buildings.HeatTransfer.Conduction.MultiLayer</a> instead of this model.
</p>
<h4>Important parameters</h4>
<p>
The parameters <code>stateAtSurface_a</code> and
<code>stateAtSurface_b</code>
determine whether there is a state variable at these surfaces,
as described above.
Note that if <code>stateAtSurface_a = true</code>,
then there is temperature state on the surface a with prescribed
value, as determined by the differential equation of the heat conduction.
Hence, in this situation, it is not possible to
connect a temperature boundary condition such as
<a href=\"modelica://Buildings.HeatTransfer.Sources.FixedTemperature\">
Buildings.HeatTransfer.Sources.FixedTemperature</a> as this would
yield to specifying the same temperature twice.
To avoid this, either set <code>stateAtSurface_a = false</code>,
or place a thermal resistance
between the boundary condition and the surface of this model.
The same applies for surface b.
See the examples in
<a href=\"modelica://Buildings.HeatTransfer.Examples\">
Buildings.HeatTransfer.Examples</a>.
</p>
</html>",
  revisions="<html>
<ul>
<li>
August 27, 2019, by Michael Wetter:<br/>
Removed assertion on geometry.<br/>
This is for
<a href=\"https://github.com/lbl-srg/modelica-buildings/issues/1529\">issue 1529</a>.
</li>
<li>
November 22, 2016, by Thierry S. Nouidui:<br/>
Fix bug in mass balance.
</li>
<li>
November 17, 2016, by Thierry S. Nouidui:<br/>
Added parameter <code>nSta2</code> to avoid translation error
in Dymola 2107. This is a work-around for a bug in Dymola
which will be addressed in future releases.
</li>
<li>
November 11, 2016, by Thierry S. Nouidui:<br/>
Revised the implementation for adding a state at the surface.
</li>
<li>
October 29, 2016, by Michael Wetter:<br/>
Added option to place a state at the surface.<br/>
This is for
<a href=\"https://github.com/lbl-srg/modelica-buildings/issues/565\">issue 565</a>.
</li>
<li>
March 1, 2016, by Michael Wetter:<br/>
Removed test for equality of <code>Real</code> variables.
This is for
<a href=\"https://github.com/lbl-srg/modelica-buildings/issues/493\">issue 493</a>.
</li>
<li>
May 21, 2015, by Michael Wetter:<br/>
Reformulated function to reduce use of the division macro
in Dymola.
This is for <a href=\"https://github.com/lbl-srg/modelica-buildings/issues/417\">issue 417</a>.
</li>
<li>
October 17, 2014, by Michael Wetter:<br/>
Changed the input argument for the function
<code>Buildings.HeatTransfer.Conduction.BaseClasses.der_temperature_u</code>
from type
<code>Buildings.HeatTransfer.Data.BaseClasses.Material</code>
to the elements of this type as OpenModelica fails to translate the
model if the input to this function is a record.
</li>
<li>
May 30, 2014, by Michael Wetter:<br/>
Removed undesirable annotation <code>Evaluate=true</code>.
</li>
<li>
January 22, 2013, by Armin Teskeredzic:<br/>
Implementation of phase-change materials based on enthalpy-linearisation method.
Phase-change properties defined in <code>material</code> record and relationship
between enthalpy and temperature defined in the <code>EnthalpyTemperature</code> function.
</li>
<li>
March 9, 2012, by Michael Wetter:<br/>
Removed protected variable <code>der_T</code> as it is not required.
</li>
<li>
March 6 2010, by Michael Wetter:<br/>
Changed implementation to allow steady-state and transient heat conduction
depending on the specific heat capacity of the material. This allows using the
same model in composite constructions in which some layers are
computed steady-state and other transient.
</li><li>
February 5 2009, by Michael Wetter:<br/>
First implementation.
</li>
</ul>
</html>"));
  end SingleLayer;

  model LectureBuilding
    "Simple thermal model of building with 4 lecture halls and a hallway"
      // Room Parameters
      parameter Modelica.SIunits.Length h_room[1,5] = [8, 8, 8, 8, 8] "Height of rooms Hall 1 -4 and Hallway";
      parameter Modelica.SIunits.Length l_room[1,5] = [30, 30, 20, 20, 50] "Length of rooms Hall 1 -4 and Hallway";
      parameter Modelica.SIunits.Length w_room[1,5] = [20, 20, 20, 20, 10] "Width of rooms Hall 1 -4 and Hallway";
      parameter Modelica.SIunits.Length h_win[1,5] = [3, 3, 3, 3, 5] "Height of Windows of Hall 1 -4 and Hallway";
      parameter Modelica.SIunits.Length w_win[1,5] = [25, 25, 15, 15, 10] "Width of Windows of Hall 1 -4 and Hallway";
      parameter Real Q_occ_rad = 75 "Radiant Heat per person [W]";
      parameter Real Q_occ_conv = 25 "Convective Heat per person [W]";
      parameter Real Q_occ_lat = 0 "Latent Heat per person [W]";

    parameter Buildings.HeatTransfer.Data.GlazingSystems.DoubleClearAir13Clear FE "Window construction"
      annotation (Placement(transformation(extent={{254,222},{274,242}})));

    Buildings.ThermalZones.Detailed.MixedAir Hs1(
      nConExt=2,
      nConExtWin=1,
      nConPar=0,
      nConBou=2,
      nSurBou=1,
      datConExt(
        layers={AW,AW},
        A={h_room[1, 1]*w_room[1, 1],w_room[1, 1]*l_room[1, 1]},
        til={Buildings.Types.Tilt.Wall,Buildings.Types.Tilt.Ceiling},
        azi={Buildings.Types.Azimuth.N,Buildings.Types.Azimuth.S}),
      datConExtWin(
        layers={AW},
        A={h_room[1, 1]*l_room[1, 1]},
        glaSys={FE},
        hWin={h_win[1, 1]},
        wWin={w_win[1, 1]},
        fFra={0.1},
        til={Buildings.Types.Tilt.Wall},
        azi={Buildings.Types.Azimuth.W}),
      datConBou(
        layers={IW,IW},
        A={l_room[1, 1]*h_room[1, 1],w_room[1, 1]*h_room[1, 1]},
        til={Buildings.Types.Tilt.Wall,Buildings.Types.Tilt.Wall}),
      surBou(
        A={l_room[1, 1]*w_room[1, 1]},
        absIR={0.9},
        absSol={0.9},
        til={Buildings.Types.Tilt.Floor}),
      redeclare package Medium = Buildings.Media.Air (T_default=293.15),
      lat=0.91664692314742,
      AFlo=l_room[1, 1]*w_room[1, 1],
      hRoo=h_room[1, 1],
      intConMod=Buildings.HeatTransfer.Types.InteriorConvection.Temperature,
      T_start=293.15,
      nPorts=2)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}},
          rotation=0,
          origin={44,100})));

    Buildings.ThermalZones.Detailed.MixedAir Flur(
      nConExt=1,
      nConExtWin=2,
      nConPar=0,
      nConBou=0,
      nSurBou=5,
      datConExt(
        layers={AW},
        A={w_room[1, 5]*l_room[1, 5]},
        til={Buildings.Types.Tilt.Ceiling},
        azi={Buildings.Types.Azimuth.S}),
      datConExtWin(
        layers={AW,AW},
        A={h_room[1, 5]*w_room[1, 5],h_room[1, 5]*w_room[1, 5]},
        glaSys={FE,FE},
        hWin={h_win[1, 5],h_win[1, 5]},
        wWin={w_win[1, 5],w_win[1, 5]},
        fFra={0.1,0.1},
        til={Buildings.Types.Tilt.Wall,Buildings.Types.Tilt.Wall},
        azi={Buildings.Types.Azimuth.N,Buildings.Types.Azimuth.S}),
      surBou(
        A={l_room[1, 1]*h_room[1, 5],l_room[1, 2]*h_room[1, 5],l_room[1, 3]*
            h_room[1, 5],l_room[1, 4]*h_room[1, 5],l_room[1, 5]*w_room[1, 5]},
        absIR={0.9,0.9,0.9,0.9,0.9},
        absSol={0.9,0.9,0.9,0.9,0.9},
        til={Buildings.Types.Tilt.Wall,Buildings.Types.Tilt.Wall,Buildings.Types.Tilt.Wall,
            Buildings.Types.Tilt.Wall,Buildings.Types.Tilt.Floor}),
      redeclare package Medium = Buildings.Media.Air (T_default=293.15),
      lat(displayUnit="deg") = 0.91664692314742,
      AFlo=l_room[1, 5]*w_room[1, 5],
      hRoo=h_room[1, 5],
      linearizeRadiation=true,
      intConMod=Buildings.HeatTransfer.Types.InteriorConvection.Temperature,
      T_start=293.15,
      nPorts=1) annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=0,
          origin={118,52})));

    Buildings.ThermalZones.Detailed.MixedAir Hs2(
      nConExt=2,
      nConExtWin=1,
      nConPar=0,
      nConBou=2,
      nSurBou=1,
      datConExt(
        layers={AW,AW},
        A={h_room[1, 2]*w_room[1, 2],w_room[1, 2]*l_room[1, 2]},
        til={Buildings.Types.Tilt.Wall,Buildings.Types.Tilt.Ceiling},
        azi={Buildings.Types.Azimuth.N,Buildings.Types.Azimuth.S}),
      datConExtWin(
        layers={AW},
        A={h_room[1, 2]*l_room[1, 2]},
        glaSys={FE},
        hWin={h_win[1, 2]},
        wWin={w_win[1, 2]},
        fFra={0.1},
        til={Buildings.Types.Tilt.Wall},
        azi={Buildings.Types.Azimuth.E}),
      datConBou(
        layers={IW,IW},
        A={l_room[1, 2]*h_room[1, 2],w_room[1, 2]*h_room[1, 2]},
        til={Buildings.Types.Tilt.Wall,Buildings.Types.Tilt.Wall}),
      surBou(
        A={l_room[1, 2]*w_room[1, 2]},
        absIR={0.9},
        absSol={0.9},
        til={Buildings.Types.Tilt.Floor}),
      redeclare package Medium = Buildings.Media.Air (T_default=293.15),
      lat=0.91664692314742,
      AFlo=l_room[1, 2]*w_room[1, 2],
      hRoo=h_room[1, 2],
      intConMod=Buildings.HeatTransfer.Types.InteriorConvection.Temperature,
      T_start=293.15,
      nPorts=2)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}},
          rotation=0,
          origin={198,100})));

    Buildings.ThermalZones.Detailed.MixedAir Hs3(
      nConExt=2,
      nConExtWin=1,
      nConPar=0,
      nConBou=1,
      nSurBou=2,
      datConExt(
        layers={AW,AW},
        A={h_room[1, 3]*w_room[1, 3],w_room[1, 3]*l_room[1, 3]},
        til={Buildings.Types.Tilt.Wall,Buildings.Types.Tilt.Ceiling},
        azi={Buildings.Types.Azimuth.S,Buildings.Types.Azimuth.S}),
      datConExtWin(
        layers={AW},
        A={h_room[1, 3]*l_room[1, 3]},
        glaSys={FE},
        hWin={h_win[1, 3]},
        wWin={w_win[1, 3]},
        fFra={0.1},
        til={Buildings.Types.Tilt.Wall},
        azi={Buildings.Types.Azimuth.W}),
      datConBou(
        layers={IW},
        A={l_room[1, 3]*h_room[1, 3]},
        til={Buildings.Types.Tilt.Wall}),
      surBou(
        A={w_room[1, 3]*h_room[1, 3],w_room[1, 3]*l_room[1, 3]},
        absIR={0.9,0.9},
        absSol={0.9,0.9},
        til={Buildings.Types.Tilt.Wall,Buildings.Types.Tilt.Floor}),
      redeclare package Medium = Buildings.Media.Air (T_default=293.15),
      lat=0.91664692314742,
      AFlo=w_room[1, 3]*l_room[1, 3],
      hRoo=h_room[1, 3],
      intConMod=Buildings.HeatTransfer.Types.InteriorConvection.Temperature,
      extConMod=Buildings.HeatTransfer.Types.ExteriorConvection.TemperatureWind,
      use_C_flow=false,
      T_start=293.15,
      nPorts=2)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}},
          rotation=0,
          origin={46,-4})));
    Buildings.ThermalZones.Detailed.MixedAir Hs4(
      nConExt=2,
      nConExtWin=1,
      nConPar=0,
      nConBou=1,
      nSurBou=2,
      datConExt(
        layers={AW,AW},
        A={h_room[1, 4]*w_room[1, 4],w_room[1, 4]*l_room[1, 4]},
        til={Buildings.Types.Tilt.Wall,Buildings.Types.Tilt.Ceiling},
        azi={Buildings.Types.Azimuth.S,Buildings.Types.Azimuth.S}),
      datConExtWin(
        layers={AW},
        A={h_room[1, 4]*l_room[1, 4]},
        glaSys={FE},
        hWin={h_win[1, 4]},
        wWin={w_win[1, 4]},
        fFra={0.1},
        til={Buildings.Types.Tilt.Wall},
        azi={Buildings.Types.Azimuth.E}),
      datConBou(
        layers={IW},
        A={l_room[1, 4]*h_room[1, 4]},
        til={Buildings.Types.Tilt.Wall}),
      surBou(
        A={w_room[1, 4]*h_room[1, 4],w_room[1, 4]*l_room[1, 4]},
        absIR={0.9,0.9},
        absSol={0.9,0.9},
        til={Buildings.Types.Tilt.Wall,Buildings.Types.Tilt.Floor}),
      redeclare package Medium = Buildings.Media.Air (T_default=293.15),
      lat=0.91664692314742,
      AFlo=w_room[1, 4]*l_room[1, 4],
      hRoo=h_room[1, 4],
      intConMod=Buildings.HeatTransfer.Types.InteriorConvection.Temperature,
      use_C_flow=false,
      T_start=293.15,
      nPorts=2)
      annotation (Placement(transformation(extent={{-20,-20},{20,20}},
          rotation=0,
          origin={200,-4})));
    Modelica.Blocks.Sources.RealExpression T0(y=-Modelica.Constants.T_zero)
                                                         "abs. Zero"
      annotation (Placement(transformation(extent={{-56,206},{-36,226}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor T_Hs1
      annotation (Placement(transformation(extent={{234,144},{254,164}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor T_Hs2
      annotation (Placement(transformation(extent={{234,118},{254,138}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor T_Hs3
      annotation (Placement(transformation(extent={{234,92},{254,112}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor T_Hs4
      annotation (Placement(transformation(extent={{234,64},{254,84}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor T_Flur annotation (
        Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={244,44})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature T_Erde
      "Temperature of Soil"
      annotation (Placement(transformation(extent={{-64,24},{-44,44}})));
    Modelica.Blocks.Interfaces.RealInput T_Sup "Supply Temperature of AHU [C]"
      annotation (Placement(transformation(extent={{-20,-20},{20,20}},
          rotation=0,
          origin={-90,204})));
    Modelica.Blocks.Sources.RealExpression T_Erde_Sig(y=9 - 6*cos(2*Modelica.Constants.pi
          *time/(86400*365) - 2*Modelica.Constants.pi*60/365) - Modelica.Constants.T_zero)
      "Temperature of Soil"
      annotation (Placement(transformation(extent={{-92,24},{-72,44}})));
    Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat(
      filNam="./DEU_Berlin.103840_IWEC.mos",
        computeWetBulbTemperature=true,
      calTSky=Buildings.BoundaryConditions.Types.SkyTemperatureCalculation.TemperaturesAndSkyCover)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=270,
          origin={136,230})));
    Buildings.Fluid.Sources.MassFlowSource_T Zul1(
      redeclare package Medium = Buildings.Media.Air,
      use_m_flow_in=true,
      m_flow=2,
      use_T_in=true,
      T=293.15,
      nPorts=1) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={16,186})));
    Buildings.Fluid.Sources.MassFlowSource_T Zul2(
      redeclare package Medium = Buildings.Media.Air,
      use_m_flow_in=true,
      m_flow=2,
      use_T_in=true,
      T=293.15,
      nPorts=1) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={40,186})));
    Buildings.Fluid.Sources.MassFlowSource_T Zul3(
      redeclare package Medium = Buildings.Media.Air,
      use_m_flow_in=true,
      m_flow=1.5,
      use_T_in=true,
      T=293.15,
      nPorts=1) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={68,186})));
    Buildings.Fluid.Sources.MassFlowSource_T Zul4(
      redeclare package Medium = Buildings.Media.Air,
      use_m_flow_in=true,
      m_flow=1.5,
      use_T_in=true,
      T=293.15,
      nPorts=1) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={94,186})));
    Modelica.Blocks.Math.Add T_ZulToK(k1=+1, k2=+1) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-10,210})));
    Modelica.Blocks.Routing.DeMultiplex4 DeMuxV(
      n1=1,
      n2=1,
      n3=1,
      n4=1) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-10,240})));
    Buildings.BoundaryConditions.WeatherData.Bus weaBus "Bus with weather data"
      annotation (Placement(transformation(extent={{126,198},{146,218}})));
    parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic FB(material={
          Buildings.HeatTransfer.Data.Solids.Concrete(x=0.3),Buildings.HeatTransfer.Data.Solids.InsulationBoard(x=0.15),
          Buildings.HeatTransfer.Data.Solids.Concrete(x=0.1)},
                                          final nLay=3) "Floor construction"
      annotation (Placement(transformation(extent={{220,222},{240,242}})));
    //Buildings.Fluid.Sources.Outside_Cp out1
      //annotation (Placement(transformation(extent={{-390,154},{-370,174}})));
    parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic IW(material={
          Buildings.HeatTransfer.Data.Solids.Concrete(x=0.2)},  final nLay=1) "Inner wall construction"
      annotation (Placement(transformation(extent={{190,222},{210,242}})));
    Modelica.Blocks.Interfaces.RealInput M_Sup[4]
      "4x1 Supply Massflowrate [kg/s]"
      annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=0,
          origin={-88,240})));
    Modelica.Blocks.Interfaces.RealInput N_Pers[4] "4x1 Number of Persons [1] "
      annotation (Placement(transformation(extent={{-20,-20},{20,20}},
          rotation=0,
          origin={-88,130})));
    Modelica.Blocks.Interfaces.RealOutput T_air[5] "5x1 Temperature [C] "
      annotation (Placement(transformation(extent={{332,106},{360,134}})));
    LecBuilding.MultiLayer BodenPlatte(
      A=Hs1.AFlo + Hs2.AFlo + Hs3.AFlo + Hs4.AFlo + Flur.AFlo,
      layers=FB,
      steadyStateInitial=false)
      annotation (Placement(transformation(extent={{-40,24},{-20,44}})));
    LecBuilding.NPersToHeatFlow nPersToHeatFlow(A=[Hs1.AFlo; Hs2.AFlo; Hs3.AFlo;
          Hs4.AFlo], Q=[Q_occ_rad; Q_occ_conv; Q_occ_lat])
      annotation (Placement(transformation(extent={{-60,120},{-40,140}})));
    Buildings.Fluid.Sources.Boundary_pT bouPres(
      redeclare package Medium = Buildings.Media.Air,
      use_p_in=false,
      p=100000,
      nPorts=5)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-76,-8})));
    Modelica.Blocks.Routing.Multiplex5 muxT annotation (Placement(transformation(extent={{306,110},
              {326,130}})));
    parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic AW(material={
          Buildings.HeatTransfer.Data.Solids.InsulationBoard(x=0.15),Buildings.HeatTransfer.Data.Solids.Concrete(x=0.3)},
        final nLay=2) "Outer wall construction" annotation (Placement(transformation(extent={{160,222},
              {180,242}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin ToC1
      annotation (Placement(transformation(extent={{264,144},{284,164}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin ToC2
      annotation (Placement(transformation(extent={{264,118},{284,138}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin ToC3
      annotation (Placement(transformation(extent={{264,92},{284,112}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin ToC4
      annotation (Placement(transformation(extent={{264,64},{284,84}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin ToC5
      annotation (Placement(transformation(extent={{264,34},{284,54}})));
    Modelica.Blocks.Math.MatrixGain gai1(K=[0; 0; 0])
      "Gain to convert from occupancy (per person) to radiant, convective and latent heat in [W/m2] "
      annotation (Placement(transformation(extent={{84,56},{92,64}})));
    Modelica.Blocks.Sources.RealExpression Zero(y=0) "Zero" annotation (
        Placement(transformation(
          extent={{-3,-9},{3,9}},
          rotation=90,
          origin={79,49})));
  equation
    connect(T_Hs2.port, Hs2.heaPorAir)
      annotation (Line(points={{234,128},{197,128},{197,100}},color={191,0,0}));
    connect(Hs3.heaPorAir, T_Hs3.port) annotation (Line(points={{45,-4},{74,-4},
            {74,-26},{226,-26},{226,102},{234,102}},
                                color={191,0,0}));
    connect(Hs4.heaPorAir, T_Hs4.port) annotation (Line(points={{199,-4},{230,
            -4},{230,74},{234,74}},
                                  color={191,0,0}));
    connect(T_Flur.port, Flur.heaPorAir) annotation (Line(points={{234,44},{176,
            44},{176,52},{117,52}},  color={191,0,0}));
    connect(T_Erde_Sig.y, T_Erde.T)
      annotation (Line(points={{-71,34},{-66,34}},            color={0,0,127}));

    connect(Zul4.ports[1], Hs4.ports[1]) annotation (Line(points={{94,176},{94,
            148},{170,148},{170,-16},{185,-16}},
                                          color={0,127,255}));
    connect(T_ZulToK.u2,T_Sup)
      annotation (Line(points={{-22,204},{-90,204}},  color={0,0,127}));
    connect(weaDat.weaBus, weaBus) annotation (Line(
        points={{136,220},{136,208}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%second",
        index=1,
        extent={{6,3},{6,3}},
        horizontalAlignment=TextAlignment.Left));
    connect(DeMuxV.u,M_Sup)  annotation (Line(points={{-22,240},{-88,240}},
                                    color={0,0,127}));
    connect(Hs2.surf_conBou[1], Flur.surf_surBou[2]) annotation (Line(points={{204,
            83.5},{204,34},{114.2,34},{114.2,37.6}},         color={191,0,0}));
    connect(Hs3.surf_conBou[1], Flur.surf_surBou[3]) annotation (Line(points={{52,-20},
            {52,24},{114.2,24},{114.2,38}},             color={191,0,0}));
    connect(Hs4.surf_conBou[1], Flur.surf_surBou[4]) annotation (Line(points={{206,-20},
            {206,24},{114.2,24},{114.2,38.4}},            color={191,0,0}));
    connect(Hs4.surf_surBou[1], Hs2.surf_conBou[2]) annotation (Line(points={{196.2,
            -18.5},{204,-18.5},{204,84.5}},                 color={191,0,0}));
    connect(Zul3.ports[1], Hs3.ports[1])
      annotation (Line(points={{68,176},{68,-26},{24,-26},{24,-16},{31,-16}},
                                                          color={0,127,255}));
    connect(BodenPlatte.port_a, T_Erde.port)
      annotation (Line(points={{-40,34},{-44,34}}, color={191,0,0}));
    connect(Hs2.surf_surBou[1], BodenPlatte.port_b) annotation (Line(points={{194.2,
            86},{194.2,34},{-20,34}},     color={191,0,0}));
    connect(Hs3.surf_surBou[2], BodenPlatte.port_b) annotation (Line(points={{42.2,
            -17.5},{50,-17.5},{50,34},{-20,34}},  color={191,0,0}));
    connect(Hs4.surf_surBou[2], BodenPlatte.port_b) annotation (Line(points={{196.2,
            -17.5},{204,-17.5},{204,34},{-20,34}},   color={191,0,0}));
    connect(Flur.surf_surBou[5], BodenPlatte.port_b) annotation (Line(points={{114.2,
            38.8},{114.2,34},{-20,34}},     color={191,0,0}));
    connect(nPersToHeatFlow.y2[:, 1], Hs2.qGai_flow) annotation (Line(points={{-39,133},
            {-40,133},{-40,132},{150,132},{150,108},{176.4,108}},
                                                             color={0,0,127}));
    connect(N_Pers, nPersToHeatFlow.u[:, 1])
      annotation (Line(points={{-88,130},{-61.8,130}},color={0,0,127}));
    connect(nPersToHeatFlow.y3[:, 1], Hs3.qGai_flow) annotation (Line(points={{-39,127},
            {-10,127},{-10,4},{24.4,4}},                    color={0,0,127}));
    connect(Hs4.qGai_flow, nPersToHeatFlow.y4[:, 1]) annotation (Line(points={{178.4,4},
            {166,4},{166,-30},{-14,-30},{-14,121},{-39,121}},  color={0,0,127}));
    connect(Zul2.ports[1], Hs2.ports[1])
      annotation (Line(points={{40,176},{40,126},{166,126},{166,88},{183,88}}, color={0,127,255}));
    connect(Hs3.ports[2], bouPres.ports[1]) annotation (Line(points={{31,-12},{
            -66,-12},{-66,-4.8}},                                                             color={0,127,255}));
    connect(DeMuxV.y1[1], Zul1.m_flow_in) annotation (Line(points={{1,249},{24,
            249},{24,198}},                                                                    color={0,0,127}));
    connect(DeMuxV.y2[1], Zul2.m_flow_in) annotation (Line(points={{1,243},{1,
            244},{48,244},{48,198}},                                                                   color={0,0,127}));
    connect(DeMuxV.y3[1], Zul3.m_flow_in) annotation (Line(points={{1,237},{76,
            237},{76,198}},                                                                    color={0,0,127}));
    connect(DeMuxV.y4[1], Zul4.m_flow_in) annotation (Line(points={{1,231},{102,
            231},{102,198}},                                                                     color={0,0,127}));
    connect(T_ZulToK.y, Zul1.T_in) annotation (Line(points={{1,210},{20,210},{
            20,198}},                                                                   color={0,0,127}));
    connect(Zul2.T_in, Zul1.T_in) annotation (Line(points={{44,198},{44,210},{
            20,210},{20,198}},                                                                   color={0,0,127}));
    connect(Zul3.T_in, Zul1.T_in) annotation (Line(points={{72,198},{72,210},{
            20,210},{20,198}},                                                                   color={0,0,127}));
    connect(Zul4.T_in, Zul1.T_in) annotation (Line(points={{98,198},{98,210},{
            20,210},{20,198}},                                                                   color={0,0,127}));
    connect(muxT.y, T_air) annotation (Line(points={{327,120},{346,120}},                     color={0,0,127}));
    connect(T_Hs1.T, ToC1.Kelvin)
      annotation (Line(points={{254,154},{262,154}}, color={0,0,127}));
    connect(ToC1.Celsius, muxT.u1[1]) annotation (Line(points={{285,154},{292,
            154},{292,130},{304,130}},
                                  color={0,0,127}));
    connect(T_Hs2.T, ToC2.Kelvin)
      annotation (Line(points={{254,128},{262,128}}, color={0,0,127}));
    connect(ToC2.Celsius, muxT.u2[1]) annotation (Line(points={{285,128},{288,
            128},{288,125},{304,125}},
                                  color={0,0,127}));
    connect(T_Hs3.T, ToC3.Kelvin) annotation (Line(points={{254,102},{262,102}},
                             color={0,0,127}));
    connect(ToC3.Celsius, muxT.u3[1]) annotation (Line(points={{285,102},{292,
            102},{292,120},{304,120}},
                                  color={0,0,127}));
    connect(T_Hs4.T, ToC4.Kelvin)
      annotation (Line(points={{254,74},{262,74}},   color={0,0,127}));
    connect(ToC4.Celsius, muxT.u4[1]) annotation (Line(points={{285,74},{296,74},
            {296,115},{304,115}}, color={0,0,127}));
    connect(T_Flur.T, ToC5.Kelvin)
      annotation (Line(points={{254,44},{262,44}}, color={0,0,127}));
    connect(ToC5.Celsius, muxT.u5[1]) annotation (Line(points={{285,44},{300,44},
            {300,110},{304,110}},color={0,0,127}));
    connect(T0.y, T_ZulToK.u1)
      annotation (Line(points={{-35,216},{-22,216}}, color={0,0,127}));
    connect(Hs1.surf_surBou[1], BodenPlatte.port_b) annotation (Line(points={{40.2,86},
            {40,86},{40,34},{-20,34}},          color={191,0,0}));
    connect(Hs1.surf_conBou[1], Flur.surf_surBou[1]) annotation (Line(points={{50,83.5},
            {50,34},{114.2,34},{114.2,37.2}},          color={191,0,0}));
    connect(Hs1.surf_conBou[2], Hs3.surf_surBou[1]) annotation (Line(points={{50,84.5},
            {50,-18.5},{42.2,-18.5}},        color={191,0,0}));
    connect(Hs1.ports[1], Zul1.ports[1])
      annotation (Line(points={{29,88},{16,88},{16,176}}, color={0,127,255}));
    connect(nPersToHeatFlow.y1[:, 1], Hs1.qGai_flow) annotation (Line(points={{-39,139},
            {-40,139},{-40,138},{-6,138},{-6,108},{22.4,108}},
                                                      color={0,0,127}));
    connect(T_Hs1.port, Hs1.heaPorAir)
      annotation (Line(points={{234,154},{43,154},{43,100}},color={191,0,0}));
    connect(Zero.y, gai1.u[1])
      annotation (Line(points={{79,52.3},{79,60},{83.2,60}}, color={0,0,127}));
    connect(gai1.y, Flur.qGai_flow)
      annotation (Line(points={{92.4,60},{96.4,60}},  color={0,0,127}));
    connect(Hs1.ports[2], bouPres.ports[2]) annotation (Line(points={{29,92},{8,
            92},{8,0},{-66,0},{-66,-6.4}}, color={0,127,255}));
    connect(Flur.ports[1], bouPres.ports[3]) annotation (Line(points={{103,42},
            {16,42},{16,-8},{-66,-8}}, color={0,127,255}));
    connect(Hs2.ports[2], bouPres.ports[4]) annotation (Line(points={{183,92},{
            162,92},{162,80},{12,80},{12,-4},{-62,-4},{-62,-6},{-66,-6},{-66,
            -9.6}}, color={0,127,255}));
    connect(weaBus, Flur.weaBus) annotation (Line(
        points={{136,208},{136,69.9},{135.9,69.9}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%first",
        index=-1,
        extent={{6,3},{6,3}},
        horizontalAlignment=TextAlignment.Left));
    connect(Hs1.weaBus, Flur.weaBus) annotation (Line(
        points={{61.9,117.9},{61.9,140},{136,140},{136,69.9},{135.9,69.9}},
        color={255,204,51},
        thickness=0.5));
    connect(Hs2.weaBus, Flur.weaBus) annotation (Line(
        points={{215.9,117.9},{215.9,128},{216,128},{216,140},{136,140},{136,
            69.9},{135.9,69.9}},
        color={255,204,51},
        thickness=0.5));
    connect(Hs4.ports[2], bouPres.ports[5]) annotation (Line(points={{185,-12},
            {178,-12},{178,-34},{-66,-34},{-66,-11.2}}, color={0,127,255}));
    connect(Hs4.weaBus, Flur.weaBus) annotation (Line(
        points={{217.9,13.9},{217.9,76},{136,76},{136,69.9},{135.9,69.9}},
        color={255,204,51},
        thickness=0.5));
    connect(Hs3.weaBus, Flur.weaBus) annotation (Line(
        points={{63.9,13.9},{63.9,76},{136,76},{136,69.9},{135.9,69.9}},
        color={255,204,51},
        thickness=0.5));
    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-40},{360,
              260}})),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-40},{
              360,260}}), graphics={
          Rectangle(
            extent={{-60,256},{112,168}},
            lineColor={28,108,200},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Rectangle(
            extent={{120,256},{300,200}},
            lineColor={28,108,200},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Rectangle(
            extent={{0,160},{230,-36}},
            lineColor={28,108,200},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Rectangle(
            extent={{-96,66},{-6,20}},
            lineColor={28,108,200},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Text(
            extent={{168,216},{294,204}},
            lineColor={0,0,0},
            pattern=LinePattern.None,
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            fontSize=18,
            fontName="Arial",
            textStyle={TextStyle.Bold},
            textString="Wetter und Materialparameter"),
          Text(
            extent={{44,8},{196,-6}},
            lineColor={0,0,0},
            pattern=LinePattern.None,
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            fontSize=18,
            fontName="Arial",
            textStyle={TextStyle.Bold},
            textString="Flur, Hörsäle, 
Decke und Wände"),
          Text(
            extent={{-60,182},{-24,172}},
            lineColor={0,0,0},
            pattern=LinePattern.None,
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            fontSize=18,
            fontName="Arial",
            textStyle={TextStyle.Bold},
            textString="Zuluft"),
          Text(
            extent={{-96,66},{-62,52}},
            lineColor={0,0,0},
            pattern=LinePattern.None,
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            fontSize=18,
            fontName="Arial",
            textStyle={TextStyle.Bold},
            textString="Boden"),
          Rectangle(
            extent={{-96,14},{-6,-36}},
            lineColor={28,108,200},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Text(
            extent={{-158,-30},{-4,-34}},
            lineColor={0,0,0},
            pattern=LinePattern.None,
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            fontSize=18,
            fontName="Arial",
            textString="Abluft",
            textStyle={TextStyle.Bold})}),
      experiment(
        StopTime=2678400,
        Tolerance=1e-06,
        __Dymola_fixedstepsize=1,
        __Dymola_Algorithm="Dassl"),
      Documentation(info="<html>
      <p>
      This example represents a simplified model of a lecture building with 4 lecture halls (Hs1 - Hs4) that are supplied 
      by the same AHU. The mass flow rates to each hall can be adjusted by using input M_Sup. Since all halls are supplied
      by the same AHU only one temperature of all massflowrates can be chosen using input Tsup. Further, using input N_Pers,     
      the number of occupants in all rooms can be adjusted (60+20+20W radiant, convective and latent heat per person). 
      Further there is a hallway that can be heated with a direct heat flow that is split up 
      into 60% radiative and 40% convective heat. Weather Data is used from Berlin. The output of the model are the air 
      temperatures of lecture hall 1 to 4 and the hallway that are collected into the output vector T_air.</p>
      <p>
      The default dimensions of the building are as follows: 
      <ul>
        <li> Lecture hall 1: 30 x 20 x 8 m^3, 2 Outside Walls facing W and N, 1 window (3 x 25) facing W, 2 inner walls
        bordering hallway (E) and lecture hall 3 (S) <\li>
        <li> Lecture hall 2: 30 x 20 x 8 m^3, 2 Outside Walls facing E and N, 1 window (3 x 25) facing E, 2 inner walls
        bordering hallway (W) and lecture hall 4 (S) <\li>
        <li> Lecture hall 3: 20 x 20 x 8 m^3, 2 Outside Walls facing W and S, 1 window (3 x 15) facing W, 2 inner walls
        bordering hallway (E) and lecture hall 1 (N) <\li>
        <li> Lecture hall 4: 20 x 20 x 8 m^3, 2 Outside Walls facing E and S, 1 window (3 x 15) facing E, 2 inner walls
        bordering hallway (W) and lecture hall 2 (N) <\li>
        <li> Hallway: 50 x 10 x 8 m^3, 2 Outside Walls facing N and S, 2 windows (both 3 x 10) facing N and S, 4 inner walls
        bordering the lecture halls <\li>
      <\ul>
      <\p>
      <p>
      All outer walls and roofs are made of 15cm insulation board and 30cm of concrete, whereas the inner walls consist of 
      20cm concrete. All thermal zones are connected to a floor (50 x 50) that is modeled outside the zones that consists of   
      30cm concrete, 15cm insulation board and 10cm conrete. The soil temperature is modeled by a sine with mean of 9, 
      amplitude of 6 degrees and frequency of 1 year (Minimum at 1st of march). 
      </p>
      </html>"));
  end LectureBuilding;
  annotation (uses(Buildings(version="7.0.0"), Modelica(version="3.2.3"),
      ModelicaServices(version="3.2.3")));
end LecBuilding;
