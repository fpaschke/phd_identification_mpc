within ;
package SingleRoom "Models a Room with floor and radiator heating"

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
    SingleRoom.SingleLayer[nLay] lay(
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
            for i in 2:(nSta-1) loop
              T[i] =T_a_start + (T_b_start - T_a_start)*UA*sum(RNod[k] for k in 1:i-1);
            end for;
          else // stateAtSurface_a == false
            for i in 1:(nSta-1) loop
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

  model TestRoom "Tests Radiator heating of the room"
    Modelica.Blocks.Sources.Pulse nPer(
      width=9/24*100,
      amplitude=10,
      nperiod=-1,
      startTime(displayUnit="h") = 32400,
      period(displayUnit="d") = 86400) "Number of persons"
      annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
    Modelica.Blocks.Continuous.LimPID PID_rad(
      controllerType=Modelica.Blocks.Types.SimpleController.PI,
      k=1,
      Ti(displayUnit="h") = 1800,
      yMax=1,
      yMin=0,
      initType=Modelica.Blocks.Types.InitPID.InitialState,
      xi_start=0)
      annotation (Placement(transformation(extent={{-46,-22},{-26,-2}})));
    Modelica.Blocks.Sources.Pulse Tset(
      width=50,
      amplitude=6,
      nperiod=10,
      offset=16,
      startTime(displayUnit="h") = 28800,
      period(displayUnit="d") = 86400) "Temperature Setpoint"
      annotation (Placement(transformation(extent={{-80,-22},{-60,-2}})));
    Modelica.Blocks.Sources.Step step(
      height=6,
      offset=16,
      startTime(displayUnit="d") = 1036800)
      annotation (Placement(transformation(extent={{-80,-62},{-60,-42}})));
    Modelica.Blocks.Continuous.LimPID PID_flo(
      controllerType=Modelica.Blocks.Types.SimpleController.PI,
      k=1,
      Ti(displayUnit="h") = 7200,
      yMax=1,
      yMin=0,
      initType=Modelica.Blocks.Types.InitPID.InitialState,
      xi_start=0)
      annotation (Placement(transformation(extent={{-20,-42},{0,-22}})));
    Modelica.Blocks.Sources.Step step1(
      height=0,
      offset=0,
      startTime(displayUnit="d") = 0)
      annotation (Placement(transformation(extent={{-80,74},{-60,94}})));
    Modelica.Blocks.Sources.IntegerStep integerStep(height=3, startTime(
          displayUnit="d") = 13374720)
      annotation (Placement(transformation(extent={{-80,10},{-60,30}})));
    Room room annotation (Placement(transformation(extent={{18,-6},{54,34}})));
  equation
    connect(Tset.y, PID_rad.u_s)
      annotation (Line(points={{-59,-12},{-48,-12}},
                                                 color={0,0,127}));
    connect(step.y, PID_flo.u_s) annotation (Line(points={{-59,-52},{-44,-52},{
            -44,-32},{-22,-32}}, color={0,0,127}));
    connect(PID_flo.u_m, PID_rad.u_m) annotation (Line(points={{-10,-44},{-10,
            -70},{-36,-70},{-36,-24}}, color={0,0,127}));
    connect(integerStep.y, room.H_cool) annotation (Line(points={{-59,20},{-22,
            20},{-22,18},{16,18}}, color={255,127,0}));
    connect(step1.y, room.H_win) annotation (Line(points={{-59,84},{-22,84},{
            -22,30},{16,30}}, color={0,0,127}));
    connect(nPer.y, room.N_Pers) annotation (Line(points={{-59,50},{-22,50},{
            -22,24},{16,24}}, color={0,0,127}));
    connect(PID_rad.y, room.H_rad) annotation (Line(points={{-25,-12},{-4,-12},
            {-4,12},{16,12}}, color={0,0,127}));
    connect(PID_flo.y, room.H_flo) annotation (Line(points={{1,-32},{8,-32},{8,
            6},{16,6}}, color={0,0,127}));
    connect(room.T_room, PID_rad.u_m) annotation (Line(points={{56.1,25.3},{80,
            25.3},{80,-70},{-36,-70},{-36,-24}}, color={0,0,127}));
    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false)),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StopTime=1728000,
        __Dymola_fixedstepsize=15,
        __Dymola_Algorithm="Dassl"));
  end TestRoom;

  model ParallelCircuitsSlab
    "Model of multiple parallel circuits of a radiant slab"
    extends Buildings.Fluid.Interfaces.PartialTwoPort(
      port_a(p(start=p_start,
               nominal=Medium.p_default)),
      port_b(p(start=p_start,
             nominal=Medium.p_default)));
    extends Buildings.Fluid.HeatExchangers.RadiantSlabs.BaseClasses.Slab;
    extends Buildings.Fluid.Interfaces.LumpedVolumeDeclarations;
    extends Buildings.Fluid.Interfaces.TwoPortFlowResistanceParameters(
     dp_nominal = Modelica.Fluid.Pipes.BaseClasses.WallFriction.Detailed.pressureLoss_m_flow(
        m_flow=m_flow_nominal/nCir,
        rho_a=rho_default,
        rho_b=rho_default,
        mu_a=mu_default,
        mu_b=mu_default,
        length=length,
        diameter=pipe.dIn,
        roughness=pipe.roughness,
        m_flow_small=m_flow_small/nCir));

    constant Boolean homotopyInitialization = true "= true, use homotopy method"
      annotation(HideResult=true);

    parameter Integer nCir(min=1) = 1 "Number of parallel circuits";
    parameter Integer nSeg(min=1) = if heatTransfer==Buildings.Fluid.HeatExchangers.RadiantSlabs.Types.HeatTransfer.EpsilonNTU
                                                                                   then 1 else 5
      "Number of volume segments in each circuit (along flow path)";

    parameter Modelica.SIunits.Area A
      "Surface area of radiant slab (all circuits combined)"
    annotation(Dialog(group="Construction"));
    parameter Modelica.SIunits.Length length = A/disPip/nCir
      "Length of the pipe of a single circuit";

    parameter Modelica.SIunits.MassFlowRate m_flow_nominal
      "Nominal mass flow rate of all circuits combined"
      annotation(Dialog(group = "Nominal condition"));
    parameter Modelica.SIunits.MassFlowRate m_flow_small(min=0) = 1E-4*abs(m_flow_nominal)
      "Small mass flow rate of all circuits combined for regularization of zero flow"
      annotation(Dialog(tab = "Advanced"));

    final parameter Modelica.SIunits.Velocity v_nominal=
      4*m_flow_nominal/pipe.dIn^2/Modelica.Constants.pi/rho_default/nCir
      "Velocity at m_flow_nominal";

    // Parameters used for the fluid model implementation

    parameter Buildings.Fluid.HeatExchangers.RadiantSlabs.Types.HeatTransfer
      heatTransfer=Buildings.Fluid.HeatExchangers.RadiantSlabs.Types.HeatTransfer.EpsilonNTU
      "Model for heat transfer between fluid and slab";

    // Diagnostics
     parameter Boolean show_T = false
      "= true, if actual temperature at port is computed"
      annotation(Dialog(tab="Advanced",group="Diagnostics"));

    Modelica.SIunits.MassFlowRate m_flow(start=0) = port_a.m_flow
      "Mass flow rate from port_a to port_b (m_flow > 0 is design flow direction) for all circuits combined";
    Modelica.SIunits.PressureDifference dp(start=0, displayUnit="Pa") = port_a.p - port_b.p
      "Pressure difference between port_a and port_b";

    Medium.ThermodynamicState sta_a=if homotopyInitialization then
        Medium.setState_phX(port_a.p,
                            homotopy(actual=noEvent(actualStream(port_a.h_outflow)),
                                     simplified=inStream(port_a.h_outflow)),
                            homotopy(actual=noEvent(actualStream(port_a.Xi_outflow)),
                                     simplified=inStream(port_a.Xi_outflow)))
      else
        Medium.setState_phX(port_a.p,
                            noEvent(actualStream(port_a.h_outflow)),
                            noEvent(actualStream(port_a.Xi_outflow))) if
           show_T "Medium properties in port_a";

    Medium.ThermodynamicState sta_b=if homotopyInitialization then
        Medium.setState_phX(port_b.p,
                            homotopy(actual=noEvent(actualStream(port_b.h_outflow)),
                                     simplified=port_b.h_outflow),
                            homotopy(actual=noEvent(actualStream(port_b.Xi_outflow)),
                              simplified=port_b.Xi_outflow))
      else
        Medium.setState_phX(port_b.p,
                            noEvent(actualStream(port_b.h_outflow)),
                            noEvent(actualStream(port_b.Xi_outflow))) if
            show_T "Medium properties in port_b";

    SingleRoom.SingleCircuitSlab sla(
      redeclare final package Medium = Medium,
      final heatTransfer=heatTransfer,
      final sysTyp=sysTyp,
      final A=A/nCir,
      final disPip=disPip,
      final pipe=pipe,
      final layers=layers,
      final steadyStateInitial=steadyStateInitial,
      final iLayPip=iLayPip,
      final T_a_start=T_a_start,
      final T_b_start=T_b_start,
      final energyDynamics=energyDynamics,
      final massDynamics=massDynamics,
      final p_start=p_start,
      final T_start=T_start,
      final X_start=X_start,
      final C_start=C_start,
      final C_nominal=C_nominal,
      final allowFlowReversal=allowFlowReversal,
      final m_flow_nominal=m_flow_nominal/nCir,
      final m_flow_small=m_flow_small/nCir,
      final homotopyInitialization=homotopyInitialization,
      final from_dp=from_dp,
      final dp_nominal=dp_nominal,
      final linearizeFlowResistance=linearizeFlowResistance,
      final deltaM=deltaM,
      final nSeg=nSeg,
      final length=length,
      final ReC=4000,
      final stateAtSurface_a=stateAtSurface_a,
      final stateAtSurface_b=stateAtSurface_b) "Single parallel circuit of the radiant slab"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
  protected
    parameter Medium.ThermodynamicState state_default = Medium.setState_pTX(
        T=Medium.T_default,
        p=Medium.p_default,
        X=Medium.X_default[1:Medium.nXi]) "Start state";
    parameter Modelica.SIunits.Density rho_default = Medium.density(state_default);
    parameter Modelica.SIunits.DynamicViscosity mu_default = Medium.dynamicViscosity(state_default)
      "Dynamic viscosity at nominal condition";

    Buildings.Fluid.BaseClasses.MassFlowRateMultiplier masFloMul_a(
        redeclare final package Medium = Medium,
        final k=nCir)
      "Mass flow multiplier, used to avoid having to instanciate multiple slab models"
      annotation (Placement(transformation(extent={{-40,-10},{-60,10}})));
    Buildings.Fluid.BaseClasses.MassFlowRateMultiplier masFloMul_b(
        redeclare final package Medium = Medium,
        final k=nCir)
      "Mass flow multiplier, used to avoid having to instanciate multiple slab models"
      annotation (Placement(transformation(extent={{40,-10},{60,10}})));
    Buildings.Fluid.HeatExchangers.RadiantSlabs.BaseClasses.HeatFlowRateMultiplier heaFloMul_a(
        final k=nCir)
      "Heat flow rate multiplier, used to avoid having to instanciate multiple slab models"
      annotation (Placement(transformation(extent={{-40,20},{-60,40}})));
    Buildings.Fluid.HeatExchangers.RadiantSlabs.BaseClasses.HeatFlowRateMultiplier heaFloMul_b(
       final k=nCir)
      "Heat flow rate multiplier, used to avoid having to instanciate multiple slab models"
      annotation (Placement(transformation(extent={{40,-40},{60,-20}})));

  initial equation
    assert(homotopyInitialization, "In " + getInstanceName() +
      ": The constant homotopyInitialization has been modified from its default value. This constant will be removed in future releases.",
      level = AssertionLevel.warning);

  equation
    connect(sla.port_b, masFloMul_b.port_a) annotation (Line(
        points={{10,0},{28,0},{28,0},{40,0}},
        color={0,127,255},
        smooth=Smooth.None));

    connect(masFloMul_b.port_b, port_b) annotation (Line(
        points={{60,0},{80,0},{80,0},{100,0}},
        color={0,127,255},
        smooth=Smooth.None));

    connect(port_a, masFloMul_a.port_b) annotation (Line(
        points={{-100,0},{-78,0},{-78,0},{-60,0}},
        color={0,127,255},
        smooth=Smooth.None));

    connect(masFloMul_a.port_a, sla.port_a) annotation (Line(
        points={{-40,0},{-24,0},{-24,0},{-10,0}},
        color={0,127,255},
        smooth=Smooth.None));

    connect(sla.surf_a,heaFloMul_a. port_a) annotation (Line(
        points={{4,10},{4,30},{-40,30}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(heaFloMul_a.port_b, surf_a) annotation (Line(
        points={{-60,30},{-70,30},{-70,50},{40,50},{40,100}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(sla.surf_b,heaFloMul_b. port_a) annotation (Line(
        points={{4,-10},{4,-30},{40,-30}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(heaFloMul_b.port_b, surf_b) annotation (Line(
        points={{60,-30},{70,-30},{70,-80},{40,-80},{40,-100}},
        color={191,0,0},
        smooth=Smooth.None));
    annotation (Documentation(info="<html>
<p>
This is a model of a radiant slab with pipes or a capillary heat exchanger
embedded in the construction.
The model is a composition of multiple models of
<a href=\"Buildings.Fluid.HeatExchangers.RadiantSlabs.SingleCircuitSlab\">
Buildings.Fluid.HeatExchangers.RadiantSlabs.SingleCircuitSlab</a>
that are arranged in a parallel.
</p>
<p>
The parameter <code>nCir</code> declares the number of parallel flow circuits.
Each circuit will have the same mass flow rate, and it is exposed to the same
port variables for the heat port at the two surfaces, and for the flow inlet and outlet.
</p>
<p>
A typical model application is as follows: Suppose a large room has a radiant slab with two parallel circuits
with the same pipe spacing and pipe length. Then, rather than using two instances of
<a href=\"Buildings.Fluid.HeatExchangers.RadiantSlabs.SingleCircuitSlab\">
Buildings.Fluid.HeatExchangers.RadiantSlabs.SingleCircuitSlab</a>,
this system can be modeled using one instance of this model in order to reduce computing effort.
See
<a href=\"modelica://Buildings.Fluid.HeatExchangers.RadiantSlabs.Examples.SingleCircuitMultipleCircuitEpsilonNTU\">
Buildings.Fluid.HeatExchangers.RadiantSlabs.Examples.SingleCircuitMultipleCircuitEpsilonNTU</a> for an example
that shows that the models give identical results.
</p>
<p>
Since this model is a parallel arrangment of <code>nCir</code> models of
<a href=\"Buildings.Fluid.HeatExchangers.RadiantSlabs.SingleCircuitSlab\">
Buildings.Fluid.HeatExchangers.RadiantSlabs.SingleCircuitSlab</a>,
we refer to
<a href=\"Buildings.Fluid.HeatExchangers.RadiantSlabs.SingleCircuitSlab\">
Buildings.Fluid.HeatExchangers.RadiantSlabs.SingleCircuitSlab</a>
for the model documentation.
</p>
<p>
See the
<a href=\"modelica://Buildings.Fluid.HeatExchangers.RadiantSlabs.UsersGuide\">
user's guide</a> for more information.
</p>
<h4>Implementation</h4>
<p>
To allow a better comment for the nominal mass flow rate, i.e., to specify that
its value is for all circuits combined, this
model does not inherit
<a href=\"modelica://Buildings.Fluid.Interfaces.PartialTwoPortInterface\">
Buildings.Fluid.Interfaces.PartialTwoPortInterface</a>.
</p>
</html>",   revisions="<html>
<ul>
<li>
April 14, 2020, by Michael Wetter:<br/>
Changed <code>homotopyInitialization</code> to a constant.<br/>
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1341\">IBPSA, #1341</a>.
</li>
<li>
January 22, 2016, by Michael Wetter:<br/>
Corrected type declaration of pressure difference.
This is
for <a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/404\">#404</a>.
</li>
<li>
June 9, 2015 by Michael Wetter:<br/>
Changed base class from
<a href=\"modelica://Modelica.Fluid.Interfaces.PartialTwoPort\">
Modelica.Fluid.Interfaces.PartialTwoPort</a>
to
<a href=\"modelica://Buildings.Fluid.Interfaces.PartialTwoPort\">
Buildings.Fluid.Interfaces.PartialTwoPort</a>.
</li>
<li>
October 10, 2013 by Michael Wetter:<br/>
Added <code>noEvent</code> to the computation of the states at the port.
This is correct, because the states are only used for reporting, but not
to compute any other variable.
Use of the states to compute other variables would violate the Modelica
language, as conditionally removed variables must not be used in any equation.
</li>
<li>
October 8, 2013, by Michael Wetter:<br/>
Removed parameter <code>show_V_flow</code>.
</li>
<li>
September 14, 2013, by Michael Wetter:<br/>
Corrected assignment of start value for pressure at <code>port_a</code>
and <code>port_b</code>, which used <code>Medium.p_default</code>
instead of the parameter <code>p_start</code>.
</li>
<li>
June 27, 2012, by Michael Wetter:<br/>
First implementation.
</li>
</ul>
</html>"),
      Icon(graphics={
          Rectangle(
            extent={{-80,80},{80,-80}},
            lineColor={95,95,95},
            lineThickness=1,
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Line(
            points={{-90,0},{-74,0},{-74,72},{60,72},{60,46},{-48,46},{-48,14},{59.8945,
                14},{59.9999,0.0136719},{92,0}},
            color={0,128,255},
            thickness=1,
            smooth=Smooth.None),
          Line(
            points={{-90,0},{-74,0},{-74,-72},{60,-72},{60,-48},{-48,-48},{-48,-12},
                {59.7891,-12},{60,0}},
            color={0,128,255},
            thickness=1,
            smooth=Smooth.None)}));
  end ParallelCircuitsSlab;

  model SingleCircuitSlab "Model of a single circuit of a radiant slab"
    extends Buildings.Fluid.HeatExchangers.RadiantSlabs.BaseClasses.Slab;
    extends Buildings.Fluid.FixedResistances.BaseClasses.Pipe(
       nSeg=if heatTransfer==Buildings.Fluid.HeatExchangers.RadiantSlabs.Types.HeatTransfer.EpsilonNTU
                                                           then 1 else 5,
       final diameter=pipe.dIn,
       length=A/disPip,
       final thicknessIns=0,
       final lambdaIns = 0.04,
       dp_nominal = Modelica.Fluid.Pipes.BaseClasses.WallFriction.Detailed.pressureLoss_m_flow(
        m_flow=m_flow_nominal,
        rho_a=rho_default,
        rho_b=rho_default,
        mu_a=mu_default,
        mu_b=mu_default,
        length=length,
        diameter=pipe.dIn,
        roughness=pipe.roughness,
        m_flow_small=m_flow_small),
        preDro(dp(nominal=200*length)));

    parameter Modelica.SIunits.Area A "Surface area of radiant slab"
      annotation(Dialog(group="Construction"));

    parameter Buildings.Fluid.HeatExchangers.RadiantSlabs.Types.HeatTransfer
      heatTransfer=Buildings.Fluid.HeatExchangers.RadiantSlabs.Types.HeatTransfer.EpsilonNTU
      "Model for heat transfer between fluid and slab";
    parameter Modelica.SIunits.Temperature T_c_start=
      (T_a_start*con_b[1].layers.R+T_b_start*con_a[1].layers.R)/layers.R
      "Initial construction temperature in the layer that contains the pipes, used if steadyStateInitial = false"
      annotation(Dialog(tab="Initialization", group="Construction"));
    final parameter Modelica.SIunits.Velocity v_nominal=
      4*m_flow_nominal/pipe.dIn^2/Modelica.Constants.pi/rho_default
      "Velocity at m_flow_nominal";

    Buildings.HeatTransfer.Conduction.MultiLayer con_a[nSeg](
      each final A=A/nSeg,
      each final steadyStateInitial=steadyStateInitial,
      each layers(
        final nLay = iLayPip,
        final material={layers.material[i] for i in 1:iLayPip},
        final absIR_a=layers.absIR_a,
        final absIR_b=layers.absIR_b,
        final absSol_a=layers.absSol_a,
        final absSol_b=layers.absSol_b,
        final roughness_a=layers.roughness_a),
      each T_a_start=T_a_start,
      each T_b_start=T_c_start,
      each stateAtSurface_a=stateAtSurface_a,
      each stateAtSurface_b=true)
      "Construction near the surface port surf_a"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=270,
          origin={40,50})));

    SingleRoom.MultiLayer con_b[nSeg](
        each final A=A/nSeg,
      each final steadyStateInitial=false,
        each layers(
          final nLay = layers.nLay-iLayPip,
          final material={layers.material[i] for i in iLayPip + 1:layers.nLay},
          final absIR_a=layers.absIR_a,
          final absIR_b=layers.absIR_b,
          final absSol_a=layers.absSol_a,
          final absSol_b=layers.absSol_b,
          final roughness_a=layers.roughness_a),
        each T_a_start=T_c_start,
        each T_b_start=T_b_start,
      each stateAtSurface_a=false,
      each stateAtSurface_b=stateAtSurface_b)
      "Construction near the surface port surf_b"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=270,
          origin={40,-58})));

  protected
    Modelica.Thermal.HeatTransfer.Components.ThermalCollector colAllToOne(
       final m=nSeg) "Connector to assign multiple heat ports to one heat port"
      annotation (Placement(transformation(
          extent={{-6,-6},{6,6}},
          origin={40,-80})));
    Modelica.Thermal.HeatTransfer.Components.ThermalCollector colAllToOne1(
       final m=nSeg) "Connector to assign multiple heat ports to one heat port"
      annotation (Placement(transformation(
          extent={{-6,-6},{6,6}},
          rotation=180,
          origin={40,76})));

    final parameter Modelica.SIunits.ThermalInsulance Rx=
        Buildings.Fluid.HeatExchangers.RadiantSlabs.BaseClasses.Functions.AverageResistance(
          disPip=disPip,
          dPipOut=pipe.dOut,
          k=layers.material[iLayPip].k,
          sysTyp=sysTyp,
          kIns=layers.material[iLayPip+1].k,
          dIns=layers.material[iLayPip+1].x)
      "Thermal insulance for average temperature in plane with pipes";

    Buildings.Fluid.HeatExchangers.RadiantSlabs.BaseClasses.PipeToSlabConductance
      fluSlaCon[nSeg](
      redeclare each final package Medium = Medium,
      each final APip=Modelica.Constants.pi*pipe.dIn*length/nSeg,
      each final RWal=Modelica.Math.log(pipe.dOut/pipe.dIn)/(2*Modelica.Constants.pi
          *pipe.k*(length/nSeg)),
      each final RFic=nSeg*Rx/A,
      each final m_flow_nominal=m_flow_nominal,
      each kc_IN_con=
          Modelica.Fluid.Dissipation.HeatTransfer.StraightPipe.kc_overall_IN_con(
          d_hyd=pipe.dIn,
          L=length/nSeg,
          K=pipe.roughness),
      each final heatTransfer=heatTransfer)
      "Conductance between fluid and the slab"
      annotation (Placement(transformation(extent={{-28,-80},{-8,-60}})));

    Modelica.SIunits.MassFraction Xi_in_a[Medium.nXi] = inStream(port_a.Xi_outflow)
      "Inflowing mass fraction at port_a";
    Modelica.SIunits.MassFraction Xi_in_b[Medium.nXi] = inStream(port_b.Xi_outflow)
      "Inflowing mass fraction at port_a";
    Modelica.Blocks.Sources.RealExpression T_a(
      final y=Medium.temperature_phX(p=port_a.p,
                                     h=inStream(port_a.h_outflow),
                                     X=cat(1,Xi_in_a,{1-sum(Xi_in_a)})))
      "Fluid temperature at port a"
      annotation (Placement(transformation(extent={{-80,-22},{-60,-2}})));
    Modelica.Blocks.Sources.RealExpression T_b(
      final y=Medium.temperature_phX(p=port_b.p,
                                     h=inStream(port_b.h_outflow),
                                     X=cat(1,Xi_in_b,{1-sum(Xi_in_b)})))
      "Fluid temperature at port b"
      annotation (Placement(transformation(extent={{-80,-36},{-60,-16}})));

    Modelica.Blocks.Sources.RealExpression mFlu_flow[nSeg](each final y=m_flow)
      "Input signal for mass flow rate"
      annotation (Placement(transformation(extent={{-80,-56},{-60,-36}})));

    Modelica.Blocks.Routing.Replicator T_a_rep(final nout=nSeg)
      "Signal replicator for T_a"
      annotation (Placement(transformation(extent={{-50,-16},{-42,-8}})));
    Modelica.Blocks.Routing.Replicator T_b_rep(final nout=nSeg)
      "Signal replicator for T_b"
      annotation (Placement(transformation(extent={{-50,-30},{-42,-22}})));
  equation
    connect(colAllToOne1.port_b,surf_a)  annotation (Line(
        points={{40,82},{40,100},{40,100}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(colAllToOne.port_b,surf_b)  annotation (Line(
        points={{40,-86},{40,-100}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(colAllToOne1.port_a, con_a.port_a) annotation (Line(
        points={{40,70},{40,60}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(colAllToOne.port_a, con_b.port_b)  annotation (Line(
        points={{40,-74},{40,-68}},
        color={191,0,0},
        smooth=Smooth.None));

    connect(fluSlaCon.fluid, vol.heatPort) annotation (Line(
        points={{-8.4,-70},{-6,-70},{-6,-28},{-1,-28}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(mFlu_flow.y, fluSlaCon.m_flow) annotation (Line(
        points={{-59,-46},{-50,-46},{-50,-66},{-29,-66}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(T_a.y, T_a_rep.u) annotation (Line(
        points={{-59,-12},{-50.8,-12}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(T_b.y, T_b_rep.u) annotation (Line(
        points={{-59,-26},{-50.8,-26}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(fluSlaCon.T_a, T_a_rep.y) annotation (Line(
        points={{-29,-60},{-36,-60},{-36,-12},{-41.6,-12}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(T_b_rep.y, fluSlaCon.T_b) annotation (Line(
        points={{-41.6,-26},{-38,-26},{-38,-63},{-29,-63}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(con_b.port_a, fluSlaCon.solid) annotation (Line(
        points={{40,-48},{40,-40},{20,-40},{20,-90},{-88,-90},{-88,-70},{-28.4,-70}},
        color={191,0,0},
        smooth=Smooth.None));

    connect(fluSlaCon.solid, con_a.port_b) annotation (Line(
        points={{-28.4,-70},{-88,-70},{-88,30},{40,30},{40,40}},
        color={191,0,0},
        smooth=Smooth.None));
    annotation (
  defaultComponentName="sla",
  Icon(graphics={
          Rectangle(
            extent={{-80,80},{80,-80}},
            lineColor={95,95,95},
            lineThickness=1,
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Line(
            points={{-90,0},{-74,0},{-74,72},{60,72},{60,34},{-62,34},{-62,-6},
                {60,-6},{60,-44},{-66,-44},{-66,-74},{74,-74}},
            color={0,128,255},
            thickness=1,
            smooth=Smooth.None),
          Line(
            points={{74,-74},{74,2},{92,2}},
            color={0,128,255},
            thickness=1,
            smooth=Smooth.None)}),
      Documentation(info="<html>
<p>
This is a model of a single flow circuit of a radiant slab with pipes or a capillary heat exchanger
embedded in the construction.
For a model with multiple parallel flow circuits, see
<a href=\"modelica://Buildings.Fluid.HeatExchangers.RadiantSlabs.ParallelCircuitsSlab\">
Buildings.Fluid.HeatExchangers.RadiantSlabs.ParallelCircuitsSlab</a>.
</p>
<p>
See the
<a href=\"modelica://Buildings.Fluid.HeatExchangers.RadiantSlabs.UsersGuide\">
user's guide</a> for more information.
</p>
</html>",
  revisions="<html>
<ul>
<li>
October 18, 2017, by Michael Wetter:<br/>
Removed state at surface b of <code>con_b</code>.
As this surface is connected to surface a of <code>con_a</code>, which
already has a state, the state can be removed, rather than relying
on the symbolic processor to remove one of these two states that are
directly coupled.
This is indeed required to avoid a warning about overdetermined initial equations.
</li>
<li>
January 06, 2016, by Thierry S. Nouidui:<br/>
Renamed parameter <code>nSta2</code> to <code>nSta</code>.
</li>
<li>
November 17, 2016, by Thierry S. Nouidui:<br/>
Added parameter <code>nSta2</code> to avoid translation error
in Dymola 2107. This is a work-around for a bug in Dymola
which will be addressed in future releases.
</li>
<li>
February 5, 2015, by Michael Wetter:<br/>
Renamed <code>res</code> to <code>preDro</code> for
<a href=\"https://github.com/lbl-srg/modelica-buildings/issues/292\">#292</a>.
</li>
<li>
September 12, 2014, by Michael Wetter:<br/>
Set start value for <code>hPip(fluid(T))</code> to avoid
a warning about conflicting start values in Dymola 2015 FD01.
</li>
<li>
February 27, 2013, by Michael Wetter:<br/>
Fixed bug in the assignment of the fictitious thermal resistance by replacing
<code>RFic[nSeg](each G=A/Rx)</code> with
<code>RFic[nSeg](each G=A/nSeg/Rx)</code>.
</li>
<li>
April 5, 2012, by Michael Wetter:<br/>
Revised implementation.
</li>
<li>
April 3, 2012, by Xiufeng Pang:<br/>
First implementation.
</li>
</ul>
</html>"));
  end SingleCircuitSlab;

  model Room "Single Room Model"
    // Room Parameters
    parameter Modelica.SIunits.Length h_room = 3;
    parameter Modelica.SIunits.Length l_room = 10;
    parameter Modelica.SIunits.Length w_room = 6;
    parameter Modelica.SIunits.Length h_win = 2;
    parameter Modelica.SIunits.Length w_win = 6;
    parameter Modelica.SIunits.Temperature T_room_nominal = 22+273.15;
    package HeatMedium = Buildings.Media.Water;
    replaceable package AirMedium = Buildings.Media.Air(T_default=293.15);

    // Radiator Heating Parameters
    parameter Modelica.SIunits.Power Q_flow_rad_nominal = 4000;
    parameter Modelica.SIunits.Temperature Tsup_rad_nom= 55+273.15;
    parameter Modelica.SIunits.Temperature Tret_rad_nom = 40+273.15;
    parameter Modelica.SIunits.MassFlowRate m_flow_nominal= Q_flow_rad_nominal/(Tsup_rad_nom-Tret_rad_nom)/HeatMedium.cp_const; // kg/s
    parameter Modelica.SIunits.PressureDifference dp_rad_nominal = 3000;
    parameter Modelica.SIunits.PressureDifference dp_valve_rad_nominal = 1500;
    parameter Modelica.SIunits.Pressure p_sink_rad = 100000;
    parameter Modelica.SIunits.Pressure p_source_rad = p_sink_rad+dp_rad_nominal+dp_valve_rad_nominal;

    // Floor Heating Parameters
    parameter Modelica.SIunits.Power Q_flow_flo_nominal = 4000;
    parameter Modelica.SIunits.Temperature Tsup_flo_nom = 40+273.15;
    parameter Modelica.SIunits.Temperature Tret_flo_nom = 30+273.15;
    parameter Modelica.SIunits.Length D_pipe = 0.1;
    parameter Integer N_circ = integer(ceil(l_room*w_room/D_pipe/100));
    parameter Modelica.SIunits.MassFlowRate m_flow_flo_nominal= Q_flow_flo_nominal/(Tsup_flo_nom-Tret_flo_nom)/HeatMedium.cp_const; // kg/s
    parameter Modelica.SIunits.PressureDifference dp_valve_flo_nominal = 100;
    parameter Modelica.SIunits.Pressure p_sink_flo = 100000;
    parameter Modelica.SIunits.Pressure p_source_flo = p_sink_flo+500+dp_valve_flo_nominal;

    // Cooling Parameters
    parameter Modelica.SIunits.MassFlowRate m_coolfan[4] = {0,1,2,3}/3600*w_room*l_room*h_room; // kg/s
    parameter Modelica.SIunits.Temperature Tsup_cool = 8+273.15;
    parameter Modelica.SIunits.MassFlowRate mwat_cool = 600/3600; //kg/s
    parameter Modelica.SIunits.Temperature Tret_cool_nom = 18+273.15;
    parameter Modelica.SIunits.Pressure p_sink_cool = 100000;

    Buildings.ThermalZones.Detailed.MixedAir Room(
      nConExt=0,
      nConExtWin=1,
      nConPar=0,
      nConBou=0,
      nSurBou=5,
      datConExtWin(
        layers={AW},
        A={h_room*l_room},
        glaSys={FE},
        hWin={h_win},
        wWin={w_win},
        fFra={0.1},
        til={Buildings.Types.Tilt.Wall},
        azi={Buildings.Types.Azimuth.W}),
      surBou(
        A={h_room*w_room,h_room*l_room,h_room*w_room,w_room*l_room,w_room*l_room},
        absIR={0.9,0.9,0.9,0.9,0.9},
        absSol={0.9,0.9,0.9,0.9,0.9},
        til={Buildings.Types.Tilt.Wall,Buildings.Types.Tilt.Wall,Buildings.Types.Tilt.Wall,
            Buildings.Types.Tilt.Floor,Buildings.Types.Tilt.Ceiling}),
      redeclare package Medium = Buildings.Media.Air (T_default=293.15),
      lat=0.91664692314742,
      AFlo=w_room*l_room,
      hRoo=h_room,
      nPorts=2)    annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=0,
          origin={14,8})));
    Buildings.BoundaryConditions.WeatherData.ReaderTMY3 WeaData(
      filNam="D:/Repository/prozessidentifikation/playground/Diss/Room/model/DEU_Berlin.103840_IWEC.mos",
      computeWetBulbTemperature=true,
      calTSky=Buildings.BoundaryConditions.Types.SkyTemperatureCalculation.TemperaturesAndSkyCover)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=270,
          origin={50,110})));

    Buildings.BoundaryConditions.WeatherData.Bus WeaBus "Bus with weather data"
      annotation (Placement(transformation(extent={{40,40},{60,60}})));
    Modelica.Blocks.Interfaces.RealInput N_Pers "Number of Persons [1] "
      annotation (Placement(transformation(extent={{-20,-20},{20,20}},
          rotation=0,
          origin={-220,60})));
    Modelica.Blocks.Math.MatrixGain NPersToHeatFlow(K=[75; 25; 0]/Room.AFlo)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-150,60})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature T_Soil
      "Temperature of Soil" annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-110,-222})));
    Modelica.Blocks.Sources.RealExpression T_Soil_Sig(y=9 - 6*cos(2*Modelica.Constants.pi
          *time/(86400*365) - 2*Modelica.Constants.pi*60/365) - Modelica.Constants.T_zero)
      "Temperature of Soil" annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-150,-222})));
    SingleRoom.MultiLayer WallNorth(
      A=w_room*h_room,
      layers=IW,
      steadyStateInitial=false) "Model of InnerWall" annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-72,96})));
    SingleRoom.MultiLayer WallEast(
      A=h_room*l_room,
      layers=IW,
      steadyStateInitial=false) "Model of InnerWall" annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={10,96})));
    SingleRoom.MultiLayer WallSouth(
      A=w_room*h_room,
      layers=IW,
      steadyStateInitial=false) "Model of InnerWall" annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-110,96})));
    SingleRoom.MultiLayer Ceiling(
      A=w_room*l_room,
      layers=CL,
      steadyStateInitial=false) "Model of InnerWall" annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-32,96})));
    Modelica.Blocks.Sources.CombiTimeTable TempTable(
      tableOnFile=true,
      tableName="tab",
      fileName="D:/Repository/prozessidentifikation/playground/Diss/Room/model/RoomTempData.txt",
      columns=2:6,
      smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
      extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
      "Temperatures of Neighboring rooms"
      annotation (Placement(transformation(extent={{-190,134},{-170,154}})));

    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature T_South
      "Temperature of Neighboring room" annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-110,120})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature T_East
      "Temperature of Neighboring room" annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={10,120})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature T_Ceil
      "Temperature of upper room" annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-32,120})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature T_North
      "Temperature of Neighboring room" annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-72,120})));
    Modelica.Blocks.Routing.DeMultiplex5 DeMultiplex5
      annotation (Placement(transformation(extent={{-150,134},{-130,154}})));
    parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic CL(material=
          {Buildings.HeatTransfer.Data.Solids.Concrete(x=0.25)}, final nLay=1)
      "Ceiling construction"
      annotation (Placement(transformation(extent={{82,142},{98,158}})));
    parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic IW(material=
          {Buildings.HeatTransfer.Data.Solids.GypsumBoard(x=0.025),
          Buildings.HeatTransfer.Data.Solids.Brick(x=0.115),
          Buildings.HeatTransfer.Data.Solids.GypsumBoard(x=0.025)}, final nLay=3)
                                                                              "Inner wall construction"
      annotation (Placement(transformation(extent={{62,142},{78,158}})));
    parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic AW(material=
          {Buildings.HeatTransfer.Data.Solids.Brick(x=0.115),
          Buildings.HeatTransfer.Data.Solids.InsulationBoard(x=0.2),
          Buildings.HeatTransfer.Data.Solids.GypsumBoard(x=0.025)},  final nLay=3)
                      "Outer wall construction" annotation (Placement(transformation(extent={{40,142},
              {56,158}})));
    Buildings.Fluid.HeatExchangers.Radiators.RadiatorEN442_2 Rad(
      redeclare package Medium = HeatMedium,
      nEle=5,
      fraRad=0.35,
      Q_flow_nominal(displayUnit="W") = Q_flow_rad_nominal,
      T_a_nominal=Tsup_rad_nom,
      T_b_nominal=Tret_rad_nom,
      TAir_nominal=T_room_nominal,
      n=1.25,
      dp_nominal=dp_rad_nominal)
      annotation (Placement(transformation(extent={{-94,-92},{-70,-68}})));
    Buildings.Fluid.Sources.Boundary_pT SouRad(
      p=p_source_rad,
      use_T_in=true,
      nPorts=1,
      redeclare package Medium = HeatMedium,
      use_p_in=false,
      T=Tsup_rad_nom)
      annotation (Placement(transformation(extent={{-140,-90},{-120,-70}})));
    Buildings.Fluid.Actuators.Valves.TwoWayEqualPercentage ValRad(
      redeclare package Medium = HeatMedium,
      m_flow_nominal=m_flow_nominal,
      CvData=Buildings.Fluid.Types.CvTypes.OpPoint,
      dpValve_nominal=dp_valve_rad_nominal,
      dpFixed_nominal=0)
      annotation (Placement(transformation(extent={{-116,-90},{-96,-70}})));
    Buildings.Fluid.Sources.Boundary_pT SinRad(
      redeclare package Medium = HeatMedium,
      p(displayUnit="bar") = p_sink_rad,
      T=Tret_rad_nom,
      nPorts=1)       "Sink"
      annotation (Placement(transformation(extent={{0,-90},{-20,-70}})));
    Modelica.Blocks.Interfaces.RealInput H_rad "Valveposition of radiator"
      annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=0,
          origin={-220,-60})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor T_RoomSens
      annotation (Placement(transformation(extent={{80,66},{100,86}})));
    Modelica.Blocks.Math.Add T_RoomToC(k1=+1, k2=-1)
      annotation (Placement(transformation(extent={{120,60},{140,80}})));
    Modelica.Blocks.Sources.RealExpression T0(y=-Modelica.Constants.T_zero)
      "abs. Zero"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=0,
          origin={90,-214})));
    Modelica.Blocks.Interfaces.RealOutput T_room "Room Temperature [C] "
      annotation (Placement(transformation(extent={{160,52},{202,94}})));
    Buildings.Controls.SetPoints.SupplyReturnTemperatureReset HcvRad(
      TSup_nominal(displayUnit="degC") = Tsup_rad_nom,
      TRet_nominal=Tret_rad_nom,
      TRoo_nominal=293.15,
      TOut_nominal(displayUnit="degC") = 258.15,
      TRoo=T_room_nominal,
      dTOutHeaBal(displayUnit="degC") = 0) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-170,-82})));
    Modelica.Blocks.Interfaces.RealOutput T_sup_rad
      "Radiator Supply Temperature [C] "
      annotation (Placement(transformation(extent={{160,-122},{202,-80}})));
    SingleRoom.ParallelCircuitsSlab FloHeat(
      redeclare package Medium = HeatMedium,
      sysTyp=Buildings.Fluid.HeatExchangers.RadiantSlabs.Types.SystemType.Floor,
      disPip=D_pipe,
      pipe=PIPE,
      layers=FLO,
      iLayPip=1,
      nCir=N_circ,
      A=w_room*l_room,
      m_flow_nominal=m_flow_flo_nominal)
      annotation (Placement(transformation(extent={{-94,-192},{-74,-172}})));
    Buildings.Fluid.Actuators.Valves.TwoWayEqualPercentage ValFloHeat(
      redeclare package Medium = HeatMedium,
      m_flow_nominal=m_flow_flo_nominal,
      CvData=Buildings.Fluid.Types.CvTypes.OpPoint,
      dpValve_nominal=dp_valve_flo_nominal,
      dpFixed_nominal=0)
      annotation (Placement(transformation(extent={{-116,-192},{-96,-172}})));
    Buildings.Fluid.Sources.Boundary_pT SinFlo(
      redeclare package Medium = HeatMedium,
      p(displayUnit="bar") = p_sink_flo,
      T=Tret_flo_nom,
      nPorts=1)      "Sink"
      annotation (Placement(transformation(extent={{0,-192},{-20,-172}})));
    Buildings.Fluid.Sources.Boundary_pT SouFlo(
      p=p_source_flo,
      use_T_in=true,
      redeclare package Medium = HeatMedium,
      use_p_in=false,
      T=Tsup_flo_nom,
      nPorts=1)
      annotation (Placement(transformation(extent={{-140,-192},{-120,-172}})));
    Buildings.Controls.SetPoints.SupplyReturnTemperatureReset HcvFlo(
      m=1.1,
      TSup_nominal(displayUnit="degC") = Tsup_flo_nom,
      TRet_nominal=Tret_flo_nom,
      TRoo_nominal=293.15,
      TOut_nominal(displayUnit="degC") = 258.15,
      TRoo=T_room_nominal,
      dTOutHeaBal(displayUnit="K") = 0)    annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-170,-184})));
    Modelica.Blocks.Interfaces.RealInput H_flo "Valveposition of floor heating"
      annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=0,
          origin={-220,-120})));
    Modelica.Blocks.Interfaces.RealOutput T_sup_flo
      "Floor Heating Supply Temperature [C] "
      annotation (Placement(transformation(extent={{160,-224},{202,-182}})));
    parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic FLO(material=
          {Buildings.HeatTransfer.Data.Solids.Generic(
          x=0.055,
          d=2000,
          k=1.4,
          c=1000),Buildings.HeatTransfer.Data.Solids.InsulationBoard(x=0.2),
          Buildings.HeatTransfer.Data.Solids.Concrete(x=0.25)},final nLay=3)
                                  "Floor Heating Concrete"
      annotation (Placement(transformation(extent={{128,142},{144,158}})));
    Modelica.Blocks.Math.Add T_SupRadToC(k1=+1, k2=-1)
      annotation (Placement(transformation(extent={{120,-116},{140,-96}})));
    Modelica.Blocks.Math.Add T_SupRadToC1(k1=+1, k2=-1)
      annotation (Placement(transformation(extent={{120,-218},{140,-198}})));
    Modelica.Blocks.Interfaces.RealInput H_win "Position of Window blinds [1] "
      annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=0,
          origin={-220,120})));
    parameter
      Buildings.HeatTransfer.Data.GlazingSystems.TripleClearAir13ClearAir13Clear
      FE annotation (Placement(transformation(extent={{106,142},{122,158}})));
    parameter Buildings.Fluid.Data.Pipes.PEX_RADTEST
                                     PIPE "Pipe material"
      annotation (Placement(transformation(extent={{20,142},{36,158}})));
    Buildings.Fluid.Sensors.MassFlowRate M_rad(redeclare package Medium =
          HeatMedium)
      annotation (Placement(transformation(extent={{-44,-90},{-24,-70}})));
    Modelica.Blocks.Interfaces.RealOutput P_rad
      "Heating power of radiator heating [kW] "
      annotation (Placement(transformation(extent={{160,-52},{202,-10}})));
    Modelica.Blocks.Math.Add dT_rad(k1=+1, k2=-1) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-62,-38})));
    Modelica.Blocks.Math.Product product annotation (Placement(transformation(extent={{-20,-40},{0,-20}})));
    Modelica.Blocks.Math.Gain c_w(k=HeatMedium.cp_const/1000)
                                         annotation (Placement(transformation(extent={{20,-40},{40,-20}})));
    Buildings.Fluid.Sensors.MassFlowRate M_flo(redeclare package Medium =
          HeatMedium)
      annotation (Placement(transformation(extent={{-46,-192},{-26,-172}})));
    Modelica.Blocks.Interfaces.RealOutput P_flo
      "Heating power of floor heating [kW] "
      annotation (Placement(transformation(extent={{160,-152},{202,-110}})));
    Modelica.Blocks.Math.Add dT_flo(k1=+1, k2=-1) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-66,-138})));
    Modelica.Blocks.Math.Product product1
                                         annotation (Placement(transformation(extent={{-20,
              -140},{0,-120}})));
    Modelica.Blocks.Math.Gain c_w1(k=HeatMedium.cp_const/1000)
                                         annotation (Placement(transformation(extent={{20,-140},
              {40,-120}})));
    Modelica.Fluid.Sensors.TemperatureTwoPort T_ret_flo(redeclare package
        Medium = HeatMedium)
      annotation (Placement(transformation(extent={{-70,-192},{-50,-172}})));
    Modelica.Fluid.Sensors.TemperatureTwoPort T_ret_rad(redeclare package
        Medium = HeatMedium)
      annotation (Placement(transformation(extent={{-66,-90},{-46,-70}})));
    Buildings.Fluid.HeatExchangers.DryCoilCounterFlow heaCoi(
      redeclare package Medium1 = HeatMedium,                redeclare package
        Medium2 = AirMedium,
      m1_flow_nominal=mwat_cool,
      m2_flow_nominal=Fan.massFlowRates[4],
      dp1_nominal=2000,
      dp2_nominal=200,
      UA_nominal=300,
      energyDynamics=Modelica.Fluid.Types.Dynamics.SteadyState,
      tau1=60,
      tau2=60,
      tau_m=120)
      annotation (Placement(transformation(extent={{-140,-42},{-120,-22}})));
    Buildings.Fluid.Movers.FlowControlled_m_flow Fan(
      redeclare package Medium = AirMedium,
      energyDynamics=Modelica.Fluid.Types.Dynamics.SteadyState,
      massDynamics=Modelica.Fluid.Types.Dynamics.SteadyState,
      m_flow_nominal=m_coolfan[4],
      inputType=Buildings.Fluid.Types.InputType.Stages,
      massFlowRates=m_coolfan)
      annotation (Placement(transformation(extent={{-190,-48},{-170,-28}})));
    Buildings.Fluid.Sources.MassFlowSource_T boundary(
      redeclare package Medium = HeatMedium,
      use_m_flow_in=false,                            m_flow=mwat_cool, T=
          Tsup_cool,
      nPorts=1)
      annotation (Placement(transformation(extent={{-166,-36},{-146,-16}})));
    Modelica.Fluid.Sensors.TemperatureTwoPort T_ret_cool(redeclare package
        Medium =
          HeatMedium)
      annotation (Placement(transformation(extent={{-118,-34},{-102,-18}})));
    Modelica.Blocks.Math.Add dT_cool(k1=+1, k2=-1)
                                                  annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-116,20})));
    Modelica.Blocks.Math.Product product2
                                         annotation (Placement(transformation(extent={{-80,30},
              {-60,50}})));
    Modelica.Blocks.Sources.RealExpression T1(y=Tsup_cool)
      "abs. Zero"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-140,4})));
    Modelica.Blocks.Sources.RealExpression T2(y=mwat_cool)
      "abs. Zero"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-106,46})));
    Modelica.Blocks.Math.Gain c_w2(k=HeatMedium.cp_const/1000)
                                         annotation (Placement(transformation(extent={{80,30},
              {100,50}})));
    Modelica.Blocks.Interfaces.RealOutput P_cool "Cooling power [kW] "
      annotation (Placement(transformation(extent={{160,18},{202,60}})));
    Buildings.Fluid.Sources.Boundary_pT SinCool(
      redeclare package Medium = HeatMedium,
      p(displayUnit="bar") = p_sink_cool,
      T=Tret_cool_nom,
      nPorts=1)       "Sink"
      annotation (Placement(transformation(extent={{-80,-36},{-100,-16}})));
    Modelica.Blocks.MathInteger.Sum sum(k={1,1}, nu=2) annotation (Placement(
          transformation(
          extent={{-6,-6},{6,6}},
          rotation=270,
          origin={-180,-14})));
    Modelica.Blocks.Sources.IntegerExpression integerExpression(y=1) annotation (
        Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-178,16})));
    Modelica.Blocks.Interfaces.IntegerInput H_cool "Fan Level of Cooling [0-3]"
      annotation (Placement(transformation(extent={{-240,-20},{-200,20}})));
  equation
    connect(WeaBus, Room.weaBus) annotation (Line(
        points={{50,50},{31.9,50},{31.9,25.9}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%first",
        index=-1,
        extent={{-3,6},{-3,6}},
        horizontalAlignment=TextAlignment.Right));
    connect(NPersToHeatFlow.y, Room.qGai_flow) annotation (Line(points={{-139,60},
            {-40,60},{-40,16},{-7.6,16}},                            color={0,0,127}));
    connect(NPersToHeatFlow.u[1], N_Pers) annotation (Line(points={{-162,60},{-220,
            60}},                                 color={0,0,127}));
    connect(T_South.port, WallSouth.port_a)
      annotation (Line(points={{-110,110},{-110,106}},
                                                     color={191,0,0}));
    connect(T_Ceil.port, Ceiling.port_a)
      annotation (Line(points={{-32,110},{-32,106}},
                                                  color={191,0,0}));
    connect(T_North.port, WallNorth.port_a)
      annotation (Line(points={{-72,110},{-72,106}},color={191,0,0}));
    connect(TempTable.y, DeMultiplex5.u)
      annotation (Line(points={{-169,144},{-152,144}}, color={0,0,127}));
    connect(WallEast.port_a,T_East. port)
      annotation (Line(points={{10,106},{10,110}},color={191,0,0}));
    connect(Room.heaPorRad,Rad. heatPortRad) annotation (Line(points={{13,4.2},{-79.6,4.2},{-79.6,-71.36}},
                                 color={191,0,0}));
    connect(Room.heaPorAir,Rad. heatPortCon)
      annotation (Line(points={{13,8},{-84.4,8},{-84.4,-71.36}},
                                                             color={191,0,0}));
    connect(ValRad.port_b, Rad.port_a)
      annotation (Line(points={{-96,-80},{-94,-80}}, color={0,127,255}));
    connect(H_rad, ValRad.y)
      annotation (Line(points={{-220,-60},{-106,-60},{-106,-68}},
                                                            color={0,0,127}));
    connect(SouRad.ports[1], ValRad.port_a)
      annotation (Line(points={{-120,-80},{-116,-80}}, color={0,127,255}));
    connect(Room.heaPorAir, T_RoomSens.port)
      annotation (Line(points={{13,8},{70,8},{70,76},{80,76}},
                                              color={191,0,0}));
    connect(T_RoomSens.T, T_RoomToC.u1)
      annotation (Line(points={{100,76},{118,76}},
                                                 color={0,0,127}));
    connect(T_RoomToC.y, T_room)
      annotation (Line(points={{141,70},{162,70},{162,73},{181,73}},
                                                   color={0,0,127}));
    connect(T_room, T_room)
      annotation (Line(points={{181,73},{181,73}},
                                                 color={0,0,127}));
    connect(WallSouth.port_b, Room.surf_surBou[1]) annotation (Line(points={{-110,86},
            {-110,80},{10.2,80},{10.2,-6.8}},
                                           color={191,0,0}));
    connect(Ceiling.port_b, Room.surf_surBou[4]) annotation (Line(points={{-32,86},
            {-32,80},{10,80},{10,42},{10.2,42},{10.2,-5.6}},
                                        color={191,0,0}));
    connect(HcvRad.TSup,SouRad. T_in) annotation (Line(points={{-159,-76},{-142,-76}},
                                          color={0,0,127}));
    connect(T_Soil_Sig.y, T_Soil.T)
      annotation (Line(points={{-139,-222},{-122,-222}}, color={0,0,127}));
    connect(ValFloHeat.port_b, FloHeat.port_a)
      annotation (Line(points={{-96,-182},{-94,-182}},
                                                     color={0,127,255}));
    connect(HcvRad.TOut, WeaBus.TDryBul) annotation (Line(points={{-182,-76},{-190,-76},{-190,-110},{50,-110},{50,
            50}},                              color={0,0,127}), Text(
        string="%second",
        index=1,
        extent={{-6,3},{-6,3}},
        horizontalAlignment=TextAlignment.Right));
    connect(SouFlo.ports[1], ValFloHeat.port_a)
      annotation (Line(points={{-120,-182},{-116,-182}},
                                                       color={0,127,255}));
    connect(HcvFlo.TSup, SouFlo.T_in)
      annotation (Line(points={{-159,-178},{-142,-178}},
                                                       color={0,0,127}));
    connect(ValFloHeat.y, H_flo) annotation (Line(points={{-106,-170},{-106,
            -120},{-220,-120}},
                   color={0,0,127}));
    connect(WeaData.weaBus, WeaBus) annotation (Line(
        points={{50,100},{50,50}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%second",
        index=1,
        extent={{-3,-6},{-3,-6}},
        horizontalAlignment=TextAlignment.Right));
    connect(T_SupRadToC.y, T_sup_rad) annotation (Line(points={{141,-106},{162,-106},{162,-101},{181,-101}},
                                  color={0,0,127}));
    connect(T_sup_flo, T_SupRadToC1.y) annotation (Line(points={{181,-203},{161,-203},{161,-208},{141,-208}},
                                    color={0,0,127}));
    connect(T_SupRadToC1.u1, SouFlo.T_in) annotation (Line(points={{118,-202},{-152,-202},{-152,-178},{-142,-178}},
                                          color={0,0,127}));
    connect(T_SupRadToC.u1, SouRad.T_in) annotation (Line(points={{118,-100},{-150,-100},{-150,-76},{-142,-76}},
                                         color={0,0,127}));
    connect(T0.y, T_SupRadToC1.u2)
      annotation (Line(points={{101,-214},{118,-214}}, color={0,0,127}));
    connect(T_SupRadToC.u2, T_SupRadToC1.u2) annotation (Line(points={{118,-112},{108,-112},{108,-214},{118,-214}},
                                                             color={0,0,127}));
    connect(T_RoomToC.u2, T_SupRadToC1.u2) annotation (Line(points={{118,64},{108,64},{108,-214},{118,-214}},
                                        color={0,0,127}));
    connect(Room.uSha[1], H_win) annotation (Line(points={{-7.6,26},{-20,26},{-20,
            74},{-130,74},{-130,120},{-220,120}},
                                color={0,0,127}));
    connect(T_Soil.port, FloHeat.surf_b) annotation (Line(points={{-100,-222},{
            -80,-222},{-80,-192}},
                               color={191,0,0}));
    connect(SinRad.ports[1], M_rad.port_b) annotation (Line(points={{-20,-80},{-24,
            -80}},                                                                        color={0,127,255}));
    connect(dT_rad.y, product.u1) annotation (Line(points={{-62,-27},{-62,-24},{-22,-24}}, color={0,0,127}));
    connect(HcvRad.TOut, HcvFlo.TOut) annotation (Line(points={{-182,-76},{-190,-76},{-190,-110},{50,-110},{50,-238},
            {-190,-238},{-190,-178},{-182,-178}}, color={0,0,127}));
    connect(dT_rad.u1, SouRad.T_in) annotation (Line(points={{-68,-50},{-68,-56},
            {-150,-56},{-150,-76},{-142,-76}}, color={0,0,127}));
    connect(product.u2, M_rad.m_flow) annotation (Line(points={{-22,-36},{-34,-36},
            {-34,-69}},                          color={0,0,127}));
    connect(product.y, c_w.u)
      annotation (Line(points={{1,-30},{18,-30}}, color={0,0,127}));
    connect(c_w.y, P_rad) annotation (Line(points={{41,-30},{112,-30},{112,-31},
            {181,-31}}, color={0,0,127}));
    connect(M_flo.port_b, SinFlo.ports[1])
      annotation (Line(points={{-26,-182},{-20,-182}}, color={0,127,255}));
    connect(product1.u2, M_flo.m_flow) annotation (Line(points={{-22,-136},{-36,
            -136},{-36,-171}}, color={0,0,127}));
    connect(product1.u1, dT_flo.y) annotation (Line(points={{-22,-124},{-66,
            -124},{-66,-127}}, color={0,0,127}));
    connect(dT_flo.u1, SouFlo.T_in) annotation (Line(points={{-72,-150},{-72,
            -154},{-152,-154},{-152,-178},{-142,-178}}, color={0,0,127}));
    connect(product1.y, c_w1.u) annotation (Line(points={{1,-130},{8,-130},{8,
            -130},{18,-130}}, color={0,0,127}));
    connect(c_w1.y, P_flo) annotation (Line(points={{41,-130},{102,-130},{102,-131},
            {181,-131}},       color={0,0,127}));
    connect(M_flo.port_a, T_ret_flo.port_b)
      annotation (Line(points={{-46,-182},{-50,-182}}, color={0,127,255}));
    connect(FloHeat.port_b, T_ret_flo.port_a)
      annotation (Line(points={{-74,-182},{-70,-182}}, color={0,127,255}));
    connect(T_ret_flo.T, dT_flo.u2)
      annotation (Line(points={{-60,-171},{-60,-150}}, color={0,0,127}));
    connect(M_rad.port_a, T_ret_rad.port_b)
      annotation (Line(points={{-44,-80},{-46,-80}}, color={0,127,255}));
    connect(Rad.port_b, T_ret_rad.port_a)
      annotation (Line(points={{-70,-80},{-66,-80}}, color={0,127,255}));
    connect(T_ret_rad.T, dT_rad.u2)
      annotation (Line(points={{-56,-69},{-56,-50}}, color={0,0,127}));
    connect(Fan.port_b, heaCoi.port_b2)
      annotation (Line(points={{-170,-38},{-140,-38}}, color={0,127,255}));
    connect(boundary.ports[1], heaCoi.port_a1)
      annotation (Line(points={{-146,-26},{-140,-26}}, color={0,127,255}));
    connect(heaCoi.port_b1, T_ret_cool.port_a)
      annotation (Line(points={{-120,-26},{-118,-26}}, color={0,127,255}));
    connect(T_ret_cool.T, dT_cool.u2)
      annotation (Line(points={{-110,-17.2},{-110,8}}, color={0,0,127}));
    connect(T1.y, dT_cool.u1)
      annotation (Line(points={{-129,4},{-122,4},{-122,8}}, color={0,0,127}));
    connect(T2.y, product2.u1)
      annotation (Line(points={{-95,46},{-82,46}}, color={0,0,127}));
    connect(dT_cool.y, product2.u2)
      annotation (Line(points={{-116,31},{-116,34},{-82,34}}, color={0,0,127}));
    connect(product2.y, c_w2.u)
      annotation (Line(points={{-59,40},{78,40}}, color={0,0,127}));
    connect(c_w2.y, P_cool) annotation (Line(points={{101,40},{142,40},{142,39},{181,
            39}}, color={0,0,127}));
    connect(Fan.port_a, Room.ports[1]) annotation (Line(points={{-190,-38},{-192,-38},
            {-192,-4},{-1,-4}}, color={0,127,255}));
    connect(heaCoi.port_a2, Room.ports[2]) annotation (Line(points={{-120,-38},{-106,
            -38},{-106,-6},{-2,-6},{-2,-4},{-1,-4},{-1,0}}, color={0,127,255}));
    connect(T_ret_cool.port_b, SinCool.ports[1])
      annotation (Line(points={{-102,-26},{-100,-26}}, color={0,127,255}));
    connect(integerExpression.y, sum.u[1]) annotation (Line(points={{-178,5},{-178,
            -8},{-177.9,-8}}, color={255,127,0}));
    connect(sum.y, Fan.stage)
      annotation (Line(points={{-180,-20.9},{-180,-26}}, color={255,127,0}));
    connect(H_cool, sum.u[2]) annotation (Line(points={{-220,0},{-182.1,0},{-182.1,
            -8}}, color={255,127,0}));
    connect(WallNorth.port_b, Room.surf_surBou[3]) annotation (Line(points={{
            -72,86},{-72,80},{10.2,80},{10.2,-6}}, color={191,0,0}));
    connect(WallEast.port_b, Room.surf_surBou[2]) annotation (Line(points={{10,
            86},{10,40},{10,-6.4},{10.2,-6.4}}, color={191,0,0}));
    connect(Room.surf_surBou[5], FloHeat.surf_a) annotation (Line(points={{10.2,
            -5.2},{10,-5.2},{10,-160},{-80,-160},{-80,-172}}, color={191,0,0}));
    connect(DeMultiplex5.y5[1], T_South.T) annotation (Line(points={{-129,136},
            {-110,136},{-110,132}}, color={0,0,127}));
    connect(DeMultiplex5.y4[1], T_North.T) annotation (Line(points={{-129,140},
            {-72,140},{-72,132}}, color={0,0,127}));
    connect(DeMultiplex5.y3[1], T_Ceil.T) annotation (Line(points={{-129,144},{
            -32,144},{-32,132}}, color={0,0,127}));
    connect(DeMultiplex5.y2[1], T_East.T) annotation (Line(points={{-129,148},{
            10,148},{10,132}}, color={0,0,127}));
    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-200,-240},{160,160}})),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-200,-240},{160,160}})),
      Documentation(info="<html>
      <p>
      This example represents a simplified room model that can be heated either by floor or radiant heating.
      The default size of the room is 6 x 10 x 3 m and the size of the window facing west is 2 x 6 m. The 
      room has one outer wall (with the) window that is aligned west (11.5cm brick, 20cm insulation and 
      2,5cm gypsum board). The inner walls are made of 2,5cm gypsum, 10cm insulation and 2,5cm gypsum board.
      The ceiling and floor are massive constructions. The floor boundary is the the soil temp. that is modeled
      as sine with mean of 9, amplitude of 6 degrees and frequency of 1 year (Minimum at 1st of march).
      All other (four) boundary temperatures are read from a table, which temp. values were recorded in a actual
      office building.   
      <p>
      The input H_rad, H_flo and H_win represent control signals of the supply valves of the radiator, floor 
      heating and window blinds, resp.. The supply temp. of the floor and radiator heating is determined by 
      two separate heating curves. Further, using input N_Pers the number of occupants in the room can be  
      adjusted (60+20+20W radiant, convective and latent heat per person). Weather Data is used from Berlin. 
      The outputs of the model are the room air (T_air), and radiator (T_sup_rad) and floor heating supply 
      temperatures (T_sup_flo) and the corresponding heating powers (P_rad and P_flo).  The room has a window 
      facing  west, which blinds can be controlled using H_win (0 opened 1 closed).
      </html>"),
      experiment(
        StartTime=12960000,
        StopTime=14688000,
        __Dymola_fixedstepsize=1,
        __Dymola_Algorithm="Euler"));
  end Room;

  model RoomFloConstBnds "Single Room Model with Floor Heating and Cooling"
    // Room Parameters
    parameter Modelica.SIunits.Length h_room = 3 "Height of room";
    parameter Modelica.SIunits.Length l_room = 10 "Length of Room";
    parameter Modelica.SIunits.Length w_room = 6 "Width of Room";
    parameter Modelica.SIunits.Length h_win = 2 "Height of Window";
    parameter Modelica.SIunits.Length w_win = 6 "Width of Window";
    parameter Modelica.SIunits.Temperature T_room_nominal = 22+273.15 "Nominal room temperature";
    package HeatMedium = Buildings.Media.Water "Heating Medium";
    replaceable package AirMedium = Buildings.Media.Air(T_default=293.15) "Air of Room";

    // Floor Heating Parameters
    parameter Modelica.SIunits.Power Q_flow_flo_nominal = 4000 "Nominal heating power of floor heating";
    parameter Modelica.SIunits.Temperature Tsup_flo_nom = 40+273.15 "Nominal supply temp. of floor heating";
    parameter Modelica.SIunits.Temperature Tret_flo_nom = 30+273.15 "Nominal return temp. of floor heating";
    parameter Modelica.SIunits.Length D_pipe = 0.1 "Diameter of pipe of floor heating";
    parameter Integer N_circ = integer(ceil(l_room*w_room/D_pipe/100)) "Number of parallel circuits of floor heating";
    parameter Modelica.SIunits.MassFlowRate m_flow_flo_nominal= Q_flow_flo_nominal/(Tsup_flo_nom-Tret_flo_nom)/HeatMedium.cp_const "Nominal mass flow rate of floor heating";
    parameter Modelica.SIunits.PressureDifference dp_valve_flo_nominal = 100 "Nom. pressure loss of floor heating";
    parameter Modelica.SIunits.Pressure p_sink_flo = 100000 "Pressure of floor heating sink";
    parameter Modelica.SIunits.Pressure p_source_flo = p_sink_flo+500+dp_valve_flo_nominal "Pressure of floor heating source";

    // Cooling Parameters
    parameter Real Tsup_flo_cool = 18 "Supply temp. of cooling during summerdays [C]";
    parameter Real T_cool_startday = 120 "Startday of cooling [day of year]";
    parameter Real T_cool_endday = 273 "Endday of cooling [day of year]";

    Buildings.ThermalZones.Detailed.MixedAir Roo(
      nConExt=0,
      nConExtWin=1,
      nConPar=0,
      nConBou=0,
      nSurBou=5,
      datConExtWin(
        layers={Ow},
        A={h_room*l_room},
        glaSys={Win},
        hWin={h_win},
        wWin={w_win},
        fFra={0.1},
        til={Buildings.Types.Tilt.Wall},
        azi={Buildings.Types.Azimuth.W}),
      surBou(
        A={h_room*w_room,h_room*l_room,h_room*w_room,w_room*l_room,w_room*l_room},
        absIR={0.9,0.9,0.9,0.9,0.9},
        absSol={0.9,0.9,0.9,0.9,0.9},
        til={Buildings.Types.Tilt.Wall,Buildings.Types.Tilt.Wall,Buildings.Types.Tilt.Wall,
            Buildings.Types.Tilt.Floor,Buildings.Types.Tilt.Ceiling}),
      redeclare package Medium = Buildings.Media.Air (T_default=293.15),
      lat=0.91664692314742,
      AFlo=w_room*l_room,
      hRoo=h_room) annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=0,
          origin={14,32})));

    Buildings.BoundaryConditions.WeatherData.ReaderTMY3 WeaData(
      filNam=
          "D:/Repository/prozessidentifikation/playground/Diss/Room/model/DEU_Berlin.mos",
      computeWetBulbTemperature=true,
      calTSky=Buildings.BoundaryConditions.Types.SkyTemperatureCalculation.TemperaturesAndSkyCover)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=270,
          origin={70,120})));

    Buildings.BoundaryConditions.WeatherData.Bus Wb "Bus with weather data"
      annotation (Placement(transformation(extent={{60,40},{80,60}})));
    Modelica.Blocks.Interfaces.RealInput N_Pers "Number of Persons [1] "
      annotation (Placement(transformation(extent={{-20,-20},{20,20}},
          rotation=0,
          origin={-230,20})));
    Modelica.Blocks.Math.MatrixGain ToHeatFlow(K=[75; 25; 0]/Roo.AFlo)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-170,20})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature Tsoil
      "Temperature of Soil" annotation (Placement(transformation(
          extent={{-8,-8},{8,8}},
          rotation=0,
          origin={-6,-90})));
    Modelica.Blocks.Sources.RealExpression TsoilSig(y=9 - 0*cos(2*Modelica.Constants.pi
          *time/(86400*365) - 2*Modelica.Constants.pi*60/365) - Modelica.Constants.T_zero)
      "Temperature of Soil" annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-42,-90})));
    SingleRoom.MultiLayer WallNorth(
      A=w_room*h_room,
      layers=Iw,
      steadyStateInitial=false) "Model of InnerWall" annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-72,88})));
    SingleRoom.MultiLayer WallEast(
      A=h_room*l_room,
      layers=Iw,
      steadyStateInitial=false) "Model of InnerWall" annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={10,88})));
    SingleRoom.MultiLayer WallSouth(
      A=w_room*h_room,
      layers=Iw,
      steadyStateInitial=false) "Model of InnerWall" annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-110,88})));
    SingleRoom.MultiLayer Ceiling(
      A=w_room*l_room,
      layers=Ceil,
      steadyStateInitial=false) "Model of InnerWall" annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-32,88})));

    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature Tsouth
      "Temperature of Neighboring room" annotation (Placement(transformation(
          extent={{-8,-8},{8,8}},
          rotation=270,
          origin={-110,112})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature Teast
      "Temperature of Neighboring room" annotation (Placement(transformation(
          extent={{-8,-8},{8,8}},
          rotation=270,
          origin={10,112})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature Tceil
      "Temperature of upper room" annotation (Placement(transformation(
          extent={{-8,-8},{8,8}},
          rotation=270,
          origin={-32,112})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature Tnorth
      "Temperature of Neighboring room" annotation (Placement(transformation(
          extent={{-8,-8},{8,8}},
          rotation=270,
          origin={-72,112})));
    parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic Ceil(material=
          {Buildings.HeatTransfer.Data.Solids.Concrete(x=0.20)}, final nLay=1)
      "Ceiling construction"
      annotation (Placement(transformation(extent={{112,98},{128,114}})));
    parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic Iw(material=
          {Buildings.HeatTransfer.Data.Solids.GypsumBoard(x=0.025),
          Buildings.HeatTransfer.Data.Solids.Brick(x=0.12),
          Buildings.HeatTransfer.Data.Solids.GypsumBoard(x=0.025)}, final nLay=3)
                                                                              "Inner wall construction"
      annotation (Placement(transformation(extent={{92,98},{108,114}})));
    parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic Ow(material=
          {Buildings.HeatTransfer.Data.Solids.Brick(x=0.115),
          Buildings.HeatTransfer.Data.Solids.InsulationBoard(x=0.2),
          Buildings.HeatTransfer.Data.Solids.GypsumBoard(x=0.025)},  final nLay=3)
                      "Outer wall construction" annotation (Placement(transformation(extent={{112,118},
              {128,134}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor TrooSen
      annotation (Placement(transformation(extent={{100,60},{120,80}})));
    Modelica.Blocks.Interfaces.RealOutput T_room "Room Temperature [C] "
      annotation (Placement(transformation(extent={{176,60},{196,80}})));
    SingleRoom.ParallelCircuitsSlab FloHeat(
      redeclare package Medium = HeatMedium,
      sysTyp=Buildings.Fluid.HeatExchangers.RadiantSlabs.Types.SystemType.Floor,
      disPip=D_pipe,
      pipe=Pip,
      layers=Flo,
      iLayPip=1,
      nCir=N_circ,
      A=w_room*l_room,
      m_flow_nominal=m_flow_flo_nominal)
      annotation (Placement(transformation(extent={{-4,-70},{16,-50}})));
    Buildings.Fluid.Actuators.Valves.TwoWayEqualPercentage Val(
      redeclare package Medium = HeatMedium,
      m_flow_nominal=m_flow_flo_nominal,
      CvData=Buildings.Fluid.Types.CvTypes.OpPoint,
      dpValve_nominal=dp_valve_flo_nominal,
      dpFixed_nominal=0)
      annotation (Placement(transformation(extent={{-30,-70},{-10,-50}})));
    Buildings.Fluid.Sources.Boundary_pT SinFlo(
      redeclare package Medium = HeatMedium,
      p(displayUnit="bar") = p_sink_flo,
      T=Tret_flo_nom,
      nPorts=1)      "Sink"
      annotation (Placement(transformation(extent={{140,-70},{120,-50}})));
    Buildings.Fluid.Sources.Boundary_pT Sou(
      p=p_source_flo,
      use_T_in=true,
      redeclare package Medium = HeatMedium,
      use_p_in=false,
      T=Tsup_flo_nom,
      nPorts=1)
      annotation (Placement(transformation(extent={{-58,-70},{-38,-50}})));
    Buildings.Controls.SetPoints.SupplyReturnTemperatureReset HcvFlo(
      m=1.1,
      TSup_nominal(displayUnit="degC") = Tsup_flo_nom,
      TRet_nominal=Tret_flo_nom,
      TRoo_nominal=293.15,
      TOut_nominal(displayUnit="degC") = 258.15,
      TRoo=T_room_nominal,
      dTOutHeaBal(displayUnit="K") = 0)    annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-120,-98})));
    Modelica.Blocks.Interfaces.RealOutput T_sup_flo
      "Floor Heating Supply Temperature [C] "
      annotation (Placement(transformation(extent={{176,-122},{196,-102}})));
    parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic Flo(material=
          {Buildings.HeatTransfer.Data.Solids.Generic(
          x=0.055,
          d=2000,
          k=1.4,
          c=1000),Buildings.HeatTransfer.Data.Solids.InsulationBoard(x=0.2),
          Buildings.HeatTransfer.Data.Solids.Concrete(x=0.25)}, final nLay=3)
      "Floor Heating Concrete"
      annotation (Placement(transformation(extent={{92,118},{108,134}})));
    Modelica.Blocks.Interfaces.RealInput H_win "Position of Window blinds [1] "
      annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=0,
          origin={-230,60})));
    parameter Buildings.Fluid.Data.Pipes.PEX_RADTEST Pip "Pipe material"
      annotation (Placement(transformation(extent={{132,118},{148,134}})));
    Buildings.Fluid.Sensors.MassFlowRate Mflo(redeclare package Medium =
          HeatMedium)
      annotation (Placement(transformation(extent={{80,-70},{100,-50}})));
    Modelica.Blocks.Interfaces.RealOutput P_flo
      "Heating power of floor heating [kW] "
      annotation (Placement(transformation(extent={{176,-10},{196,10}})));
    Modelica.Blocks.Math.Add DeltaTflo(k1=+1, k2=-1) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={40,-10})));
    Modelica.Blocks.Math.Product Prod
      annotation (Placement(transformation(extent={{100,-10},{120,10}})));
    Modelica.Blocks.Math.Gain Cw(k=HeatMedium.cp_const/1000)
      annotation (Placement(transformation(extent={{140,-10},{160,10}})));
    Modelica.Fluid.Sensors.TemperatureTwoPort Tret(redeclare package Medium =
          HeatMedium)
      annotation (Placement(transformation(extent={{36,-70},{56,-50}})));
    Modelica.Blocks.Logical.Switch Sw
      annotation (Placement(transformation(extent={{-118,-78},{-98,-58}})));
    Modelica.Blocks.Logical.GreaterEqualThreshold Ge(threshold=T_cool_startday*86400)
      annotation (Placement(transformation(extent={{-200,-70},{-180,-50}})));
    Modelica.Blocks.Sources.RealExpression Tcool(y=Tsup_flo_cool - Modelica.Constants.T_zero)
      "Supply Temp. Cooling" annotation (Placement(transformation(
          extent={{-9,-10},{9,10}},
          rotation=0,
          origin={-151,-60})));
    Modelica.Blocks.Logical.LessEqualThreshold Le(threshold=T_cool_endday*86400)
      annotation (Placement(transformation(extent={{-200,-110},{-180,-90}})));
    Modelica.Blocks.Logical.And And
      annotation (Placement(transformation(extent={{-160,-90},{-140,-70}})));
    Modelica.Blocks.Interfaces.RealInput H_flo "Valveposition of floor heating"
      annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=0,
          origin={-230,-20})));
    parameter Buildings.HeatTransfer.Data.GlazingSystems.DoubleClearAir13Clear Win(
        haveExteriorShade=true, shade=Buildings.HeatTransfer.Data.Shades.Generic())
      annotation (Placement(transformation(extent={{132,98},{148,114}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin ToC1
      annotation (Placement(transformation(extent={{140,60},{160,80}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin ToC2
      annotation (Placement(transformation(extent={{130,-122},{150,-102}})));
    Modelica.Blocks.Sources.RealExpression Twalls(y=20 - Modelica.Constants.T_zero)
      "Temperature of Soil" annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-168,132})));
  equation
    connect(ToHeatFlow.y, Roo.qGai_flow) annotation (Line(points={{-159,20},{
            -120,20},{-120,40},{-7.6,40}}, color={0,0,127}));
    connect(ToHeatFlow.u[1], N_Pers)
      annotation (Line(points={{-182,20},{-230,20}}, color={0,0,127}));
    connect(Tsouth.port, WallSouth.port_a)
      annotation (Line(points={{-110,104},{-110,98}},  color={191,0,0}));
    connect(Tceil.port, Ceiling.port_a)
      annotation (Line(points={{-32,104},{-32,98}},  color={191,0,0}));
    connect(Tnorth.port, WallNorth.port_a)
      annotation (Line(points={{-72,104},{-72,98}},  color={191,0,0}));
    connect(WallEast.port_a, Teast.port)
      annotation (Line(points={{10,98},{10,104}},  color={191,0,0}));
    connect(Roo.heaPorAir, TrooSen.port) annotation (Line(points={{13,32},{50,
            32},{50,70},{100,70}},
                               color={191,0,0}));
    connect(WallSouth.port_b, Roo.surf_surBou[1]) annotation (Line(points={{-110,78},
            {-110,74},{10.2,74},{10.2,17.2}}, color={191,0,0}));
    connect(Ceiling.port_b, Roo.surf_surBou[4]) annotation (Line(points={{-32,78},
            {-32,74},{10,74},{10,36},{10.2,36},{10.2,18.4}}, color={191,0,0}));
    connect(TsoilSig.y, Tsoil.T)
      annotation (Line(points={{-31,-90},{-15.6,-90}}, color={0,0,127}));
    connect(Val.port_b, FloHeat.port_a)
      annotation (Line(points={{-10,-60},{-4,-60}}, color={0,127,255}));
    connect(Sou.ports[1], Val.port_a)
      annotation (Line(points={{-38,-60},{-30,-60}}, color={0,127,255}));
    connect(WeaData.weaBus, Wb) annotation (Line(
        points={{70,110},{70,50}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%second",
        index=1,
        extent={{-3,-6},{-3,-6}},
        horizontalAlignment=TextAlignment.Right));
    connect(Roo.uSha[1], H_win) annotation (Line(points={{-7.6,50},{-120,50},{
            -120,60},{-230,60}},
                            color={0,0,127}));
    connect(Tsoil.port, FloHeat.surf_b)
      annotation (Line(points={{2,-90},{10,-90},{10,-70}},   color={191,0,0}));
    connect(Mflo.port_b, SinFlo.ports[1])
      annotation (Line(points={{100,-60},{120,-60}}, color={0,127,255}));
    connect(Prod.u2, Mflo.m_flow)
      annotation (Line(points={{98,-6},{90,-6},{90,-49}}, color={0,0,127}));
    connect(Prod.u1, DeltaTflo.y)
      annotation (Line(points={{98,6},{40,6},{40,1}}, color={0,0,127}));
    connect(Prod.y, Cw.u)
      annotation (Line(points={{121,0},{138,0}}, color={0,0,127}));
    connect(Mflo.port_a, Tret.port_b)
      annotation (Line(points={{80,-60},{56,-60}}, color={0,127,255}));
    connect(FloHeat.port_b, Tret.port_a)
      annotation (Line(points={{16,-60},{36,-60}}, color={0,127,255}));
    connect(Tret.T, DeltaTflo.u2)
      annotation (Line(points={{46,-49},{46,-22}}, color={0,0,127}));
    connect(WallNorth.port_b, Roo.surf_surBou[3]) annotation (Line(points={{-72,78},
            {-72,74},{10.2,74},{10.2,18}}, color={191,0,0}));
    connect(WallEast.port_b, Roo.surf_surBou[2])
      annotation (Line(points={{10,78},{10,17.6},{10.2,17.6}}, color={191,0,0}));
    connect(Roo.surf_surBou[5], FloHeat.surf_a) annotation (Line(points={{10.2,18.8},
            {10,18.8},{10,-50}}, color={191,0,0}));
    connect(Sw.y, Sou.T_in) annotation (Line(points={{-97,-68},{-68,-68},{-68,
            -56},{-60,-56}},
                        color={0,0,127}));
    connect(HcvFlo.TSup, Sw.u3) annotation (Line(points={{-126,-87},{-126,-76},
            {-120,-76}},
                   color={0,0,127}));
    connect(Tcool.y, Sw.u1) annotation (Line(points={{-141.1,-60},{-120,-60}},
                              color={0,0,127}));
    connect(P_flo, Cw.y)
      annotation (Line(points={{186,0},{161,0}}, color={0,0,127}));
    connect(Le.u, Wb.solTim) annotation (Line(points={{-202,-100},{-220,-100},{
            -220,-120},{70,-120},{70,50}},
                                      color={0,0,127}), Text(
        string="%second",
        index=1,
        extent={{-6,3},{-6,3}},
        horizontalAlignment=TextAlignment.Right));
    connect(Le.u, Ge.u) annotation (Line(points={{-202,-100},{-220,-100},{-220,
            -60},{-202,-60}},
                         color={0,0,127}));
    connect(Ge.y, And.u1) annotation (Line(points={{-179,-60},{-170,-60},{-170,
            -80},{-162,-80}},
                         color={255,0,255}));
    connect(Le.y, And.u2)
      annotation (Line(points={{-179,-100},{-170,-100},{-170,-88},{-162,-88}},
                                                       color={255,0,255}));
    connect(And.y, Sw.u2)
      annotation (Line(points={{-139,-80},{-130,-80},{-130,-68},{-120,-68}},
                                                       color={255,0,255}));
    connect(Val.y, H_flo) annotation (Line(points={{-20,-48},{-20,-20},{-230,
            -20}},
          color={0,0,127}));
    connect(TrooSen.T, ToC1.Kelvin)
      annotation (Line(points={{120,70},{138,70}}, color={0,0,127}));
    connect(ToC2.Kelvin, Sou.T_in) annotation (Line(points={{128,-112},{-68,
            -112},{-68,-56},{-60,-56}},
                                  color={0,0,127}));
    connect(DeltaTflo.u1, Sou.T_in) annotation (Line(points={{34,-22},{34,-40},
            {-68,-40},{-68,-56},{-60,-56}},
                                       color={0,0,127}));
    connect(Roo.weaBus, Wb) annotation (Line(
        points={{31.9,49.9},{45.95,49.9},{45.95,50},{70,50}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%second",
        index=1,
        extent={{6,3},{6,3}},
        horizontalAlignment=TextAlignment.Left));
    connect(ToC1.Celsius, T_room)
      annotation (Line(points={{161,70},{186,70}}, color={0,0,127}));
    connect(ToC2.Celsius, T_sup_flo)
      annotation (Line(points={{151,-112},{186,-112}}, color={0,0,127}));
    connect(HcvFlo.TOut, Wb.TDryBul) annotation (Line(points={{-126,-110},{-126,
            -120},{70,-120},{70,50}}, color={0,0,127}), Text(
        string="%second",
        index=1,
        extent={{-3,-6},{-3,-6}},
        horizontalAlignment=TextAlignment.Right));
    connect(Twalls.y, Tsouth.T) annotation (Line(points={{-157,132},{-110,132},
            {-110,121.6}}, color={0,0,127}));
    connect(Tnorth.T, Tsouth.T) annotation (Line(points={{-72,121.6},{-72,132},
            {-110,132},{-110,121.6}}, color={0,0,127}));
    connect(Tceil.T, Tsouth.T) annotation (Line(points={{-32,121.6},{-32,132},{
            -110,132},{-110,121.6}}, color={0,0,127}));
    connect(Teast.T, Tsouth.T) annotation (Line(points={{10,121.6},{10,132},{
            -110,132},{-110,121.6}}, color={0,0,127}));
    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-240,-140},{200,
              160}})),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-240,-140},{
              200,160}}),
                      graphics={
          Rectangle(
            extent={{-190,150},{30,68}},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Text(
            extent={{-44,146},{58,142}},
            lineColor={0,0,0},
            fillColor={255,170,170},
            fillPattern=FillPattern.None,
            fontSize=14,
            fontName="Arial",
            horizontalAlignment=TextAlignment.Left,
            textStyle={TextStyle.Bold},
            textString="Innenwände und Decke"),
          Rectangle(
            extent={{-224,-36},{-90,-130}},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Text(
            extent={{-176,-42},{-64,-46}},
            lineColor={0,0,0},
            fillColor={255,170,170},
            fillPattern=FillPattern.None,
            fontSize=14,
            fontName="Arial",
            horizontalAlignment=TextAlignment.Left,
            textStyle={TextStyle.Bold},
            textString="Vorlauftemperaturregelung"),
          Rectangle(
            extent={{-80,-36},{160,-130}},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Text(
            extent={{-74,-124},{56,-128}},
            lineColor={0,0,0},
            fillColor={255,170,170},
            fillPattern=FillPattern.None,
            fontSize=14,
            fontName="Arial",
            horizontalAlignment=TextAlignment.Left,
            textStyle={TextStyle.Bold},
            textString="Fußboden und Fußbodenheizung"),
          Rectangle(
            extent={{-80,60},{40,12}},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Text(
            extent={{-76,38},{40,6}},
            lineColor={0,0,0},
            fillColor={255,170,170},
            fillPattern=FillPattern.None,
            fontSize=14,
            fontName="Arial",
            horizontalAlignment=TextAlignment.Left,
            textStyle={TextStyle.Bold},
            textString="Raum mit Außen-
wand und Fenster"),
          Rectangle(
            extent={{44,150},{160,90}},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Text(
            extent={{52,146},{168,142}},
            lineColor={0,0,0},
            fillColor={255,170,170},
            fillPattern=FillPattern.None,
            fontSize=14,
            fontName="Arial",
            horizontalAlignment=TextAlignment.Left,
            textStyle={TextStyle.Bold},
            textString="Wetter und Materialparameter")}),
      Documentation(info="<html>
      <p>
      This example represents a simplified room model that can be heated and cooled by floor heating.
      The default size of the room is 6 x 10 x 3 m and the size of the window facing west is 2 x 6 m. The 
      room has one outer wall (with the) window that is aligned west (11.5cm brick, 20cm insulation and 
      2,5cm gypsum board). The inner walls are made of 2,5cm gypsum, 12cm brick and 2,5cm gypsum board.
      The ceiling and floor are massive constructions. The floor boundary is the the soil temp. that is modeled
      as sine with mean of 9 derees, amplitude of 6 K and frequency of 1 year (Minimum at 1st of march).
      All other (four) boundary temperatures are read from a table, which temp. values were recorded in a actual
      office building.   
      <p>
      The input H_flo and H_win represent control signals of the floor heating supply valves and window blinds, 
      resp.. The supply temp. of the floor heating is determined by heating curve in winter and is constant during 
      summer cooling. Further, using input N_Pers the number of occupants in the room can be  
      adjusted (75+25+0W radiant, convective and latent heat per person). Weather Data is used from Chemnitz. 
      The outputs of the model are the room air (T_air), and floor heating supply 
      temperatures (T_sup_flo) and the corresponding heating powers (P_flo).  The room has a window 
      facing  west, which blinds can be controlled using H_win (0 opened 1 closed).
      </html>"),
      experiment(
        StopTime=2592000,
        __Dymola_fixedstepsize=15,
        __Dymola_Algorithm="Dassl"));
  end RoomFloConstBnds;

  model RoomFloVisual "Single Room Model with Floor Heating and Cooling"
    // Room Parameters
    parameter Modelica.SIunits.Length h_room = 3 "Height of room";
    parameter Modelica.SIunits.Length l_room = 10 "Length of Room";
    parameter Modelica.SIunits.Length w_room = 6 "Width of Room";
    parameter Modelica.SIunits.Length h_win = 2 "Height of Window";
    parameter Modelica.SIunits.Length w_win = 6 "Width of Window";
    parameter Modelica.SIunits.Temperature T_room_nominal = 22+273.15 "Nominal room temperature";
    package HeatMedium = Buildings.Media.Water "Heating Medium";
    replaceable package AirMedium = Buildings.Media.Air(T_default=293.15) "Air of Room";

    // Floor Heating Parameters
    parameter Modelica.SIunits.Power Q_flow_flo_nominal = 4000 "Nominal heating power of floor heating";
    parameter Modelica.SIunits.Temperature Tsup_flo_nom = 40+273.15 "Nominal supply temp. of floor heating";
    parameter Modelica.SIunits.Temperature Tret_flo_nom = 30+273.15 "Nominal return temp. of floor heating";
    parameter Modelica.SIunits.Length D_pipe = 0.1 "Diameter of pipe of floor heating";
    parameter Integer N_circ = integer(ceil(l_room*w_room/D_pipe/100)) "Number of parallel circuits of floor heating";
    parameter Modelica.SIunits.MassFlowRate m_flow_flo_nominal= Q_flow_flo_nominal/(Tsup_flo_nom-Tret_flo_nom)/HeatMedium.cp_const "Nominal mass flow rate of floor heating";
    parameter Modelica.SIunits.PressureDifference dp_valve_flo_nominal = 100 "Nom. pressure loss of floor heating";
    parameter Modelica.SIunits.Pressure p_sink_flo = 100000 "Pressure of floor heating sink";
    parameter Modelica.SIunits.Pressure p_source_flo = p_sink_flo+500+dp_valve_flo_nominal "Pressure of floor heating source";

    // Cooling Parameters
    parameter Real Tsup_flo_cool = 18 "Supply temp. of cooling during summerdays [C]";
    parameter Real T_cool_startday = 120 "Startday of cooling [day of year]";
    parameter Real T_cool_endday = 273 "Endday of cooling [day of year]";

    Buildings.ThermalZones.Detailed.MixedAir Roo(
      nConExt=0,
      nConExtWin=1,
      nConPar=0,
      nConBou=0,
      nSurBou=5,
      datConExtWin(
        layers={Ow},
        A={h_room*l_room},
        glaSys={Win},
        hWin={h_win},
        wWin={w_win},
        fFra={0.1},
        til={Buildings.Types.Tilt.Wall},
        azi={Buildings.Types.Azimuth.W}),
      surBou(
        A={h_room*w_room,h_room*l_room,h_room*w_room,w_room*l_room,w_room*l_room},
        absIR={0.9,0.9,0.9,0.9,0.9},
        absSol={0.9,0.9,0.9,0.9,0.9},
        til={Buildings.Types.Tilt.Wall,Buildings.Types.Tilt.Wall,Buildings.Types.Tilt.Wall,
            Buildings.Types.Tilt.Floor,Buildings.Types.Tilt.Ceiling}),
      redeclare package Medium = Buildings.Media.Air (T_default=293.15),
      lat=0.91664692314742,
      AFlo=w_room*l_room,
      hRoo=h_room) annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=0,
          origin={14,32})));

    Buildings.BoundaryConditions.WeatherData.ReaderTMY3 WeaData(
      filNam=
          "D:/Repository/prozessidentifikation/playground/Diss/Room/model/DEU_Berlin.mos",
      computeWetBulbTemperature=true,
      calTSky=Buildings.BoundaryConditions.Types.SkyTemperatureCalculation.TemperaturesAndSkyCover)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=270,
          origin={70,120})));

    Buildings.BoundaryConditions.WeatherData.Bus Wb "Bus with weather data"
      annotation (Placement(transformation(extent={{60,40},{80,60}})));
    Modelica.Blocks.Interfaces.RealInput N_Pers "Number of Persons [1] "
      annotation (Placement(transformation(extent={{-20,-20},{20,20}},
          rotation=0,
          origin={-230,20})));
    Modelica.Blocks.Math.MatrixGain ToHeatFlow(K=[75; 25; 0]/Roo.AFlo)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-170,20})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature Tsoil
      "Temperature of Soil" annotation (Placement(transformation(
          extent={{-8,-8},{8,8}},
          rotation=0,
          origin={-6,-72})));
    Modelica.Blocks.Sources.RealExpression TsoilSig(y=9 - 6*cos(2*Modelica.Constants.pi
          *time/(86400*365) - 2*Modelica.Constants.pi*60/365) - Modelica.Constants.T_zero)
      "Temperature of Soil" annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-42,-72})));
    SingleRoom.MultiLayer WallNorth(
      A=w_room*h_room,
      layers=Iw,
      steadyStateInitial=false) "Model of InnerWall" annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-72,88})));
    SingleRoom.MultiLayer WallEast(
      A=h_room*l_room,
      layers=Iw,
      steadyStateInitial=false) "Model of InnerWall" annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={10,88})));
    SingleRoom.MultiLayer WallSouth(
      A=w_room*h_room,
      layers=Iw,
      steadyStateInitial=false) "Model of InnerWall" annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-110,88})));
    SingleRoom.MultiLayer Ceiling(
      A=w_room*l_room,
      layers=Ceil,
      steadyStateInitial=false) "Model of InnerWall" annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-32,88})));
    Modelica.Blocks.Sources.CombiTimeTable TemTab(
      tableOnFile=true,
      tableName="tab",
      fileName="D:/Repository/prozessidentifikation/playground/Diss/Room/model/RoomTempData.txt",
      columns=2:6,
      smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
      extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
      "Temperatures of Neighboring rooms"
      annotation (Placement(transformation(extent={{-180,124},{-160,144}})));

    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature Tsouth
      "Temperature of Neighboring room" annotation (Placement(transformation(
          extent={{-8,-8},{8,8}},
          rotation=270,
          origin={-110,112})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature Teast
      "Temperature of Neighboring room" annotation (Placement(transformation(
          extent={{-8,-8},{8,8}},
          rotation=270,
          origin={10,112})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature Tceil
      "Temperature of upper room" annotation (Placement(transformation(
          extent={{-8,-8},{8,8}},
          rotation=270,
          origin={-32,112})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature Tnorth
      "Temperature of Neighboring room" annotation (Placement(transformation(
          extent={{-8,-8},{8,8}},
          rotation=270,
          origin={-72,112})));
    Modelica.Blocks.Routing.DeMultiplex5 DeMux
      annotation (Placement(transformation(extent={{-150,124},{-130,144}})));
    parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic Ceil(material=
          {Buildings.HeatTransfer.Data.Solids.Concrete(x=0.20)}, final nLay=1)
      "Ceiling construction"
      annotation (Placement(transformation(extent={{112,98},{128,114}})));
    parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic Iw(material=
          {Buildings.HeatTransfer.Data.Solids.GypsumBoard(x=0.025),
          Buildings.HeatTransfer.Data.Solids.Brick(x=0.12),
          Buildings.HeatTransfer.Data.Solids.GypsumBoard(x=0.025)}, final nLay=3)
                                                                              "Inner wall construction"
      annotation (Placement(transformation(extent={{92,98},{108,114}})));
    parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic Ow(material=
          {Buildings.HeatTransfer.Data.Solids.Brick(x=0.115),
          Buildings.HeatTransfer.Data.Solids.InsulationBoard(x=0.2),
          Buildings.HeatTransfer.Data.Solids.GypsumBoard(x=0.025)},  final nLay=3)
                      "Outer wall construction" annotation (Placement(transformation(extent={{112,118},
              {128,134}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor TrooSen
      annotation (Placement(transformation(extent={{100,22},{120,42}})));
    Modelica.Blocks.Interfaces.RealOutput T_room "Room Temperature [C] "
      annotation (Placement(transformation(extent={{176,22},{196,42}})));
    SingleRoom.ParallelCircuitsSlab FloHeat(
      redeclare package Medium = HeatMedium,
      sysTyp=Buildings.Fluid.HeatExchangers.RadiantSlabs.Types.SystemType.Floor,
      disPip=D_pipe,
      pipe=Pip,
      layers=Flo,
      iLayPip=1,
      nCir=N_circ,
      A=w_room*l_room,
      m_flow_nominal=m_flow_flo_nominal)
      annotation (Placement(transformation(extent={{-4,-52},{16,-32}})));
    Buildings.Fluid.Actuators.Valves.TwoWayEqualPercentage Val(
      redeclare package Medium = HeatMedium,
      m_flow_nominal=m_flow_flo_nominal,
      CvData=Buildings.Fluid.Types.CvTypes.OpPoint,
      dpValve_nominal=dp_valve_flo_nominal,
      dpFixed_nominal=0)
      annotation (Placement(transformation(extent={{-30,-52},{-10,-32}})));
    Buildings.Fluid.Sources.Boundary_pT SinFlo(
      redeclare package Medium = HeatMedium,
      p(displayUnit="bar") = p_sink_flo,
      T=Tret_flo_nom,
      nPorts=1)      "Sink"
      annotation (Placement(transformation(extent={{140,-52},{120,-32}})));
    Buildings.Fluid.Sources.Boundary_pT Sou(
      p=p_source_flo,
      use_T_in=true,
      redeclare package Medium = HeatMedium,
      use_p_in=false,
      T=Tsup_flo_nom,
      nPorts=1)
      annotation (Placement(transformation(extent={{-58,-52},{-38,-32}})));
    Buildings.Controls.SetPoints.SupplyReturnTemperatureReset HcvFlo(
      m=1.1,
      TSup_nominal(displayUnit="degC") = Tsup_flo_nom,
      TRet_nominal=Tret_flo_nom,
      TRoo_nominal=293.15,
      TOut_nominal(displayUnit="degC") = 258.15,
      TRoo=T_room_nominal,
      dTOutHeaBal(displayUnit="K") = 0)    annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-120,-68})));
    parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic Flo(material=
          {Buildings.HeatTransfer.Data.Solids.Generic(
          x=0.055,
          d=2000,
          k=1.4,
          c=1000),Buildings.HeatTransfer.Data.Solids.InsulationBoard(x=0.2),
          Buildings.HeatTransfer.Data.Solids.Concrete(x=0.25)}, final nLay=3)
      "Floor Heating Concrete"
      annotation (Placement(transformation(extent={{92,118},{108,134}})));
    Modelica.Blocks.Interfaces.RealInput H_win "Position of Window blinds [1] "
      annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=0,
          origin={-230,50})));
    parameter Buildings.Fluid.Data.Pipes.PEX_RADTEST Pip "Pipe material"
      annotation (Placement(transformation(extent={{132,118},{148,134}})));
    Buildings.Fluid.Sensors.MassFlowRate Mflo(redeclare package Medium =
          HeatMedium)
      annotation (Placement(transformation(extent={{80,-52},{100,-32}})));
    Modelica.Fluid.Sensors.TemperatureTwoPort Tret(redeclare package Medium =
          HeatMedium)
      annotation (Placement(transformation(extent={{36,-52},{56,-32}})));
    Modelica.Blocks.Logical.Switch Sw
      annotation (Placement(transformation(extent={{-118,-48},{-98,-28}})));
    Modelica.Blocks.Logical.GreaterEqualThreshold Ge(threshold=T_cool_startday*86400)
      annotation (Placement(transformation(extent={{-200,-46},{-180,-26}})));
    Modelica.Blocks.Sources.RealExpression Tcool(y=Tsup_flo_cool - Modelica.Constants.T_zero)
      "Supply Temp. Cooling" annotation (Placement(transformation(
          extent={{-9,-10},{9,10}},
          rotation=0,
          origin={-151,-30})));
    Modelica.Blocks.Logical.LessEqualThreshold Le(threshold=T_cool_endday*86400)
      annotation (Placement(transformation(extent={{-200,-80},{-180,-60}})));
    Modelica.Blocks.Logical.And And
      annotation (Placement(transformation(extent={{-160,-60},{-140,-40}})));
    Modelica.Blocks.Interfaces.RealInput H_flo "Valveposition of floor heating"
      annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=0,
          origin={-230,-10})));
    parameter Buildings.HeatTransfer.Data.GlazingSystems.DoubleClearAir13Clear Win(
        haveExteriorShade=true, shade=Buildings.HeatTransfer.Data.Shades.Generic())
      annotation (Placement(transformation(extent={{132,98},{148,114}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin ToC1
      annotation (Placement(transformation(extent={{140,22},{160,42}})));
  equation
    connect(ToHeatFlow.y, Roo.qGai_flow) annotation (Line(points={{-159,20},{
            -120,20},{-120,40},{-7.6,40}}, color={0,0,127}));
    connect(ToHeatFlow.u[1], N_Pers)
      annotation (Line(points={{-182,20},{-230,20}}, color={0,0,127}));
    connect(Tsouth.port, WallSouth.port_a)
      annotation (Line(points={{-110,104},{-110,98}},  color={191,0,0}));
    connect(Tceil.port, Ceiling.port_a)
      annotation (Line(points={{-32,104},{-32,98}},  color={191,0,0}));
    connect(Tnorth.port, WallNorth.port_a)
      annotation (Line(points={{-72,104},{-72,98}},  color={191,0,0}));
    connect(TemTab.y, DeMux.u)
      annotation (Line(points={{-159,134},{-152,134}}, color={0,0,127}));
    connect(WallEast.port_a, Teast.port)
      annotation (Line(points={{10,98},{10,104}},  color={191,0,0}));
    connect(Roo.heaPorAir, TrooSen.port) annotation (Line(points={{13,32},{100,
            32}},              color={191,0,0}));
    connect(WallSouth.port_b, Roo.surf_surBou[1]) annotation (Line(points={{-110,78},
            {-110,74},{10.2,74},{10.2,17.2}}, color={191,0,0}));
    connect(Ceiling.port_b, Roo.surf_surBou[4]) annotation (Line(points={{-32,78},
            {-32,74},{10,74},{10,36},{10.2,36},{10.2,18.4}}, color={191,0,0}));
    connect(TsoilSig.y, Tsoil.T)
      annotation (Line(points={{-31,-72},{-15.6,-72}}, color={0,0,127}));
    connect(Val.port_b, FloHeat.port_a)
      annotation (Line(points={{-10,-42},{-4,-42}}, color={0,127,255}));
    connect(Sou.ports[1], Val.port_a)
      annotation (Line(points={{-38,-42},{-30,-42}}, color={0,127,255}));
    connect(WeaData.weaBus, Wb) annotation (Line(
        points={{70,110},{70,50}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%second",
        index=1,
        extent={{-3,-6},{-3,-6}},
        horizontalAlignment=TextAlignment.Right));
    connect(Roo.uSha[1], H_win) annotation (Line(points={{-7.6,50},{-230,50}},
                            color={0,0,127}));
    connect(Tsoil.port, FloHeat.surf_b)
      annotation (Line(points={{2,-72},{10,-72},{10,-52}},   color={191,0,0}));
    connect(Mflo.port_b, SinFlo.ports[1])
      annotation (Line(points={{100,-42},{120,-42}}, color={0,127,255}));
    connect(Mflo.port_a, Tret.port_b)
      annotation (Line(points={{80,-42},{56,-42}}, color={0,127,255}));
    connect(FloHeat.port_b, Tret.port_a)
      annotation (Line(points={{16,-42},{36,-42}}, color={0,127,255}));
    connect(WallNorth.port_b, Roo.surf_surBou[3]) annotation (Line(points={{-72,78},
            {-72,74},{10.2,74},{10.2,18}}, color={191,0,0}));
    connect(WallEast.port_b, Roo.surf_surBou[2])
      annotation (Line(points={{10,78},{10,17.6},{10.2,17.6}}, color={191,0,0}));
    connect(Roo.surf_surBou[5], FloHeat.surf_a) annotation (Line(points={{10.2,
            18.8},{10,18.8},{10,-32}},
                                 color={191,0,0}));
    connect(Sw.y, Sou.T_in) annotation (Line(points={{-97,-38},{-60,-38}},
                        color={0,0,127}));
    connect(HcvFlo.TSup, Sw.u3) annotation (Line(points={{-126,-57},{-126,-46},
            {-120,-46}},
                   color={0,0,127}));
    connect(Tcool.y, Sw.u1) annotation (Line(points={{-141.1,-30},{-120,-30}},
                              color={0,0,127}));
    connect(Le.u, Wb.solTim) annotation (Line(points={{-202,-70},{-220,-70},{
            -220,-90},{70,-90},{70,50}},
                                      color={0,0,127}), Text(
        string="%second",
        index=1,
        extent={{-6,3},{-6,3}},
        horizontalAlignment=TextAlignment.Right));
    connect(Le.u, Ge.u) annotation (Line(points={{-202,-70},{-220,-70},{-220,
            -36},{-202,-36}},
                         color={0,0,127}));
    connect(Ge.y, And.u1) annotation (Line(points={{-179,-36},{-170,-36},{-170,
            -50},{-162,-50}},
                         color={255,0,255}));
    connect(Le.y, And.u2)
      annotation (Line(points={{-179,-70},{-170,-70},{-170,-58},{-162,-58}},
                                                       color={255,0,255}));
    connect(And.y, Sw.u2)
      annotation (Line(points={{-139,-50},{-130,-50},{-130,-38},{-120,-38}},
                                                       color={255,0,255}));
    connect(Val.y, H_flo) annotation (Line(points={{-20,-30},{-20,-10},{-230,
            -10}},
          color={0,0,127}));
    connect(DeMux.y5[1], Tsouth.T) annotation (Line(points={{-129,126},{-110,
            126},{-110,121.6}},
                         color={0,0,127}));
    connect(DeMux.y4[1], Tnorth.T) annotation (Line(points={{-129,130},{-72,130},
            {-72,121.6}},
                       color={0,0,127}));
    connect(DeMux.y3[1], Tceil.T) annotation (Line(points={{-129,134},{-32,134},
            {-32,121.6}},
                   color={0,0,127}));
    connect(DeMux.y2[1], Teast.T)
      annotation (Line(points={{-129,138},{10,138},{10,121.6}},
                                                              color={0,0,127}));
    connect(TrooSen.T, ToC1.Kelvin)
      annotation (Line(points={{120,32},{138,32}}, color={0,0,127}));
    connect(Roo.weaBus, Wb) annotation (Line(
        points={{31.9,49.9},{45.95,49.9},{45.95,50},{70,50}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%second",
        index=1,
        extent={{6,3},{6,3}},
        horizontalAlignment=TextAlignment.Left));
    connect(ToC1.Celsius, T_room)
      annotation (Line(points={{161,32},{186,32}}, color={0,0,127}));
    connect(HcvFlo.TOut, Wb.TDryBul) annotation (Line(points={{-126,-80},{-126,
            -90},{70,-90},{70,50}}, color={0,0,127}), Text(
        string="%second",
        index=1,
        extent={{-3,-6},{-3,-6}},
        horizontalAlignment=TextAlignment.Right));
    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-240,-120},{200,
              160}})),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-240,-120},{
              200,160}}),
                      graphics={
          Rectangle(
            extent={{-190,150},{30,68}},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Text(
            extent={{-44,146},{58,142}},
            lineColor={0,0,0},
            fillColor={255,170,170},
            fillPattern=FillPattern.None,
            fontSize=14,
            fontName="Arial",
            horizontalAlignment=TextAlignment.Left,
            textStyle={TextStyle.Bold},
            textString="Innenwände und Decke"),
          Rectangle(
            extent={{-226,-20},{-90,-100}},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Text(
            extent={{-170,-94},{-58,-98}},
            lineColor={0,0,0},
            fillColor={255,170,170},
            fillPattern=FillPattern.None,
            fontSize=14,
            fontName="Arial",
            horizontalAlignment=TextAlignment.Left,
            textStyle={TextStyle.Bold},
            textString="Vorlauftemperaturregelung"),
          Rectangle(
            extent={{-80,-20},{160,-100}},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Text(
            extent={{-74,-94},{56,-98}},
            lineColor={0,0,0},
            fillColor={255,170,170},
            fillPattern=FillPattern.None,
            fontSize=14,
            fontName="Arial",
            horizontalAlignment=TextAlignment.Left,
            textStyle={TextStyle.Bold},
            textString="Fußboden und Fußbodenheizung"),
          Rectangle(
            extent={{-80,60},{40,12}},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Text(
            extent={{-76,38},{40,6}},
            lineColor={0,0,0},
            fillColor={255,170,170},
            fillPattern=FillPattern.None,
            fontSize=14,
            fontName="Arial",
            horizontalAlignment=TextAlignment.Left,
            textStyle={TextStyle.Bold},
            textString="Raum mit Außen-
wand und Fenster"),
          Rectangle(
            extent={{44,150},{160,90}},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Text(
            extent={{52,146},{168,142}},
            lineColor={0,0,0},
            fillColor={255,170,170},
            fillPattern=FillPattern.None,
            fontSize=14,
            fontName="Arial",
            horizontalAlignment=TextAlignment.Left,
            textStyle={TextStyle.Bold},
            textString="Wetter und Materialparameter")}),
      Documentation(info="<html>
      <p>
      This example represents a simplified room model that can be heated and cooled by floor heating.
      The default size of the room is 6 x 10 x 3 m and the size of the window facing west is 2 x 6 m. The 
      room has one outer wall (with the) window that is aligned west (11.5cm brick, 20cm insulation and 
      2,5cm gypsum board). The inner walls are made of 2,5cm gypsum, 12cm brick and 2,5cm gypsum board.
      The ceiling and floor are massive constructions. The floor boundary is the the soil temp. that is modeled
      as sine with mean of 9 derees, amplitude of 6 K and frequency of 1 year (Minimum at 1st of march).
      All other (four) boundary temperatures are read from a table, which temp. values were recorded in a actual
      office building.   
      <p>
      The input H_flo and H_win represent control signals of the floor heating supply valves and window blinds, 
      resp.. The supply temp. of the floor heating is determined by heating curve in winter and is constant during 
      summer cooling. Further, using input N_Pers the number of occupants in the room can be  
      adjusted (75+25+0W radiant, convective and latent heat per person). Weather Data is used from Chemnitz. 
      The outputs of the model are the room air (T_air), and floor heating supply 
      temperatures (T_sup_flo) and the corresponding heating powers (P_flo).  The room has a window 
      facing  west, which blinds can be controlled using H_win (0 opened 1 closed).
      </html>"),
      experiment(
        StopTime=2592000,
        __Dymola_fixedstepsize=15,
        __Dymola_Algorithm="Dassl"));
  end RoomFloVisual;

  model RoomFlo "Single Room Model with Floor Heating and Cooling"
    // Room Parameters
    parameter Modelica.SIunits.Length h_room = 3 "Height of room";
    parameter Modelica.SIunits.Length l_room = 10 "Length of Room";
    parameter Modelica.SIunits.Length w_room = 6 "Width of Room";
    parameter Modelica.SIunits.Length h_win = 2 "Height of Window";
    parameter Modelica.SIunits.Length w_win = 6 "Width of Window";
    parameter Modelica.SIunits.Temperature T_room_nominal = 22+273.15 "Nominal room temperature";
    package HeatMedium = Buildings.Media.Water "Heating Medium";
    replaceable package AirMedium = Buildings.Media.Air(T_default=293.15) "Air of Room";

    // Floor Heating Parameters
    parameter Modelica.SIunits.Power Q_flow_flo_nominal = 4000 "Nominal heating power of floor heating";
    parameter Modelica.SIunits.Temperature Tsup_flo_nom = 40+273.15 "Nominal supply temp. of floor heating";
    parameter Modelica.SIunits.Temperature Tret_flo_nom = 30+273.15 "Nominal return temp. of floor heating";
    parameter Modelica.SIunits.Length D_pipe = 0.1 "Diameter of pipe of floor heating";
    parameter Integer N_circ = integer(ceil(l_room*w_room/D_pipe/100)) "Number of parallel circuits of floor heating";
    parameter Modelica.SIunits.MassFlowRate m_flow_flo_nominal= Q_flow_flo_nominal/(Tsup_flo_nom-Tret_flo_nom)/HeatMedium.cp_const "Nominal mass flow rate of floor heating";
    parameter Modelica.SIunits.PressureDifference dp_valve_flo_nominal = 100 "Nom. pressure loss of floor heating";
    parameter Modelica.SIunits.Pressure p_sink_flo = 100000 "Pressure of floor heating sink";
    parameter Modelica.SIunits.Pressure p_source_flo = p_sink_flo+500+dp_valve_flo_nominal "Pressure of floor heating source";

    // Cooling Parameters
    parameter Real Tsup_flo_cool = 18 "Supply temp. of cooling during summerdays [C]";
    parameter Real T_cool_startday = 120 "Startday of cooling [day of year]";
    parameter Real T_cool_endday = 273 "Endday of cooling [day of year]";

    Buildings.ThermalZones.Detailed.MixedAir Roo(
      nConExt=0,
      nConExtWin=1,
      nConPar=0,
      nConBou=0,
      nSurBou=5,
      datConExtWin(
        layers={Ow},
        A={h_room*l_room},
        glaSys={Win},
        hWin={h_win},
        wWin={w_win},
        fFra={0.1},
        til={Buildings.Types.Tilt.Wall},
        azi={Buildings.Types.Azimuth.W}),
      surBou(
        A={h_room*w_room,h_room*l_room,h_room*w_room,w_room*l_room,w_room*l_room},
        absIR={0.9,0.9,0.9,0.9,0.9},
        absSol={0.9,0.9,0.9,0.9,0.9},
        til={Buildings.Types.Tilt.Wall,Buildings.Types.Tilt.Wall,Buildings.Types.Tilt.Wall,
            Buildings.Types.Tilt.Floor,Buildings.Types.Tilt.Ceiling}),
      redeclare package Medium = Buildings.Media.Air (T_default=293.15),
      lat=0.91664692314742,
      AFlo=w_room*l_room,
      hRoo=h_room) annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=0,
          origin={14,32})));

    Buildings.BoundaryConditions.WeatherData.ReaderTMY3 WeaData(
      filNam="./DEU_Berlin.103840_IWEC.mos",
      computeWetBulbTemperature=true,
      calTSky=Buildings.BoundaryConditions.Types.SkyTemperatureCalculation.TemperaturesAndSkyCover)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=270,
          origin={70,120})));

    Buildings.BoundaryConditions.WeatherData.Bus Wb "Bus with weather data"
      annotation (Placement(transformation(extent={{60,40},{80,60}})));
    Modelica.Blocks.Interfaces.RealInput N_Pers "Number of Persons [1] "
      annotation (Placement(transformation(extent={{-20,-20},{20,20}},
          rotation=0,
          origin={-230,20})));
    Modelica.Blocks.Math.MatrixGain ToHeatFlow(K=[75; 25; 0]/Roo.AFlo)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-170,20})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature Tsoil
      "Temperature of Soil" annotation (Placement(transformation(
          extent={{-8,-8},{8,8}},
          rotation=0,
          origin={-6,-90})));
    Modelica.Blocks.Sources.RealExpression TsoilSig(y=9 - 6*cos(2*Modelica.Constants.pi
          *time/(86400*365) - 2*Modelica.Constants.pi*60/365) - Modelica.Constants.T_zero)
      "Temperature of Soil" annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-42,-90})));
    SingleRoom.MultiLayer WallNorth(
      A=w_room*h_room,
      layers=Iw,
      steadyStateInitial=false) "Model of InnerWall" annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-72,88})));
    SingleRoom.MultiLayer WallEast(
      A=h_room*l_room,
      layers=Iw,
      steadyStateInitial=false) "Model of InnerWall" annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={10,88})));
    SingleRoom.MultiLayer WallSouth(
      A=w_room*h_room,
      layers=Iw,
      steadyStateInitial=false) "Model of InnerWall" annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-110,88})));
    SingleRoom.MultiLayer Ceiling(
      A=w_room*l_room,
      layers=Ceil,
      steadyStateInitial=false) "Model of InnerWall" annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-32,88})));
    Modelica.Blocks.Sources.CombiTimeTable TemTab(
      tableOnFile=true,
      tableName="tab",
      fileName="./RoomTempData.txt",
      columns=2:6,
      smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
      extrapolation=Modelica.Blocks.Types.Extrapolation.Periodic)
      "Temperatures of Neighboring rooms"
      annotation (Placement(transformation(extent={{-180,124},{-160,144}})));

    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature Tsouth
      "Temperature of Neighboring room" annotation (Placement(transformation(
          extent={{-8,-8},{8,8}},
          rotation=270,
          origin={-110,112})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature Teast
      "Temperature of Neighboring room" annotation (Placement(transformation(
          extent={{-8,-8},{8,8}},
          rotation=270,
          origin={10,112})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature Tceil
      "Temperature of upper room" annotation (Placement(transformation(
          extent={{-8,-8},{8,8}},
          rotation=270,
          origin={-32,112})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature Tnorth
      "Temperature of Neighboring room" annotation (Placement(transformation(
          extent={{-8,-8},{8,8}},
          rotation=270,
          origin={-72,112})));
    Modelica.Blocks.Routing.DeMultiplex5 DeMux
      annotation (Placement(transformation(extent={{-150,124},{-130,144}})));
    parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic Ceil(material=
          {Buildings.HeatTransfer.Data.Solids.Concrete(x=0.20)}, final nLay=1)
      "Ceiling construction"
      annotation (Placement(transformation(extent={{112,98},{128,114}})));
    parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic Iw(material=
          {Buildings.HeatTransfer.Data.Solids.GypsumBoard(x=0.025),
          Buildings.HeatTransfer.Data.Solids.Brick(x=0.12),
          Buildings.HeatTransfer.Data.Solids.GypsumBoard(x=0.025)}, final nLay=3)
                                                                              "Inner wall construction"
      annotation (Placement(transformation(extent={{92,98},{108,114}})));
    parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic Ow(material=
          {Buildings.HeatTransfer.Data.Solids.Brick(x=0.115),
          Buildings.HeatTransfer.Data.Solids.InsulationBoard(x=0.2),
          Buildings.HeatTransfer.Data.Solids.GypsumBoard(x=0.025)},  final nLay=3)
                      "Outer wall construction" annotation (Placement(transformation(extent={{112,118},
              {128,134}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor TrooSen
      annotation (Placement(transformation(extent={{100,60},{120,80}})));
    Modelica.Blocks.Interfaces.RealOutput T_room "Room Temperature [C] "
      annotation (Placement(transformation(extent={{176,60},{196,80}})));
    SingleRoom.ParallelCircuitsSlab FloHeat(
      redeclare package Medium = HeatMedium,
      sysTyp=Buildings.Fluid.HeatExchangers.RadiantSlabs.Types.SystemType.Floor,
      disPip=D_pipe,
      pipe=Pip,
      layers=Flo,
      iLayPip=1,
      nCir=N_circ,
      A=w_room*l_room,
      m_flow_nominal=m_flow_flo_nominal)
      annotation (Placement(transformation(extent={{-4,-70},{16,-50}})));
    Buildings.Fluid.Actuators.Valves.TwoWayEqualPercentage Val(
      redeclare package Medium = HeatMedium,
      m_flow_nominal=m_flow_flo_nominal,
      CvData=Buildings.Fluid.Types.CvTypes.OpPoint,
      dpValve_nominal=dp_valve_flo_nominal,
      dpFixed_nominal=0)
      annotation (Placement(transformation(extent={{-30,-70},{-10,-50}})));
    Buildings.Fluid.Sources.Boundary_pT SinFlo(
      redeclare package Medium = HeatMedium,
      p(displayUnit="bar") = p_sink_flo,
      T=Tret_flo_nom,
      nPorts=1)      "Sink"
      annotation (Placement(transformation(extent={{140,-70},{120,-50}})));
    Buildings.Fluid.Sources.Boundary_pT Sou(
      p=p_source_flo,
      use_T_in=true,
      redeclare package Medium = HeatMedium,
      use_p_in=false,
      T=Tsup_flo_nom,
      nPorts=1)
      annotation (Placement(transformation(extent={{-58,-70},{-38,-50}})));
    Buildings.Controls.SetPoints.SupplyReturnTemperatureReset HcvFlo(
      m=1.1,
      TSup_nominal(displayUnit="degC") = Tsup_flo_nom,
      TRet_nominal=Tret_flo_nom,
      TRoo_nominal=293.15,
      TOut_nominal(displayUnit="degC") = 258.15,
      TRoo=T_room_nominal,
      dTOutHeaBal(displayUnit="K") = 0)    annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-120,-98})));
    Modelica.Blocks.Interfaces.RealOutput T_sup_flo
      "Floor Heating Supply Temperature [C] "
      annotation (Placement(transformation(extent={{176,-122},{196,-102}})));
    parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic Flo(material=
          {Buildings.HeatTransfer.Data.Solids.Generic(
          x=0.055,
          d=2000,
          k=1.4,
          c=1000),Buildings.HeatTransfer.Data.Solids.InsulationBoard(x=0.2),
          Buildings.HeatTransfer.Data.Solids.Concrete(x=0.25)}, final nLay=3)
      "Floor Heating Concrete"
      annotation (Placement(transformation(extent={{92,118},{108,134}})));
    Modelica.Blocks.Interfaces.RealInput H_win "Position of Window blinds [1] "
      annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=0,
          origin={-230,60})));
    parameter Buildings.Fluid.Data.Pipes.PEX_RADTEST Pip "Pipe material"
      annotation (Placement(transformation(extent={{132,118},{148,134}})));
    Buildings.Fluid.Sensors.MassFlowRate Mflo(redeclare package Medium =
          HeatMedium)
      annotation (Placement(transformation(extent={{80,-70},{100,-50}})));
    Modelica.Blocks.Interfaces.RealOutput P_flo
      "Heating power of floor heating [kW] "
      annotation (Placement(transformation(extent={{176,-10},{196,10}})));
    Modelica.Blocks.Math.Add DeltaTflo(k1=+1, k2=-1) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={40,-10})));
    Modelica.Blocks.Math.Product Prod
      annotation (Placement(transformation(extent={{100,-10},{120,10}})));
    Modelica.Blocks.Math.Gain Cw(k=HeatMedium.cp_const/1000)
      annotation (Placement(transformation(extent={{140,-10},{160,10}})));
    Modelica.Fluid.Sensors.TemperatureTwoPort Tret(redeclare package Medium =
          HeatMedium)
      annotation (Placement(transformation(extent={{36,-70},{56,-50}})));
    Modelica.Blocks.Logical.Switch Sw
      annotation (Placement(transformation(extent={{-118,-78},{-98,-58}})));
    Modelica.Blocks.Logical.GreaterEqualThreshold Ge(threshold=T_cool_startday*86400)
      annotation (Placement(transformation(extent={{-200,-70},{-180,-50}})));
    Modelica.Blocks.Sources.RealExpression Tcool(y=Tsup_flo_cool - Modelica.Constants.T_zero)
      "Supply Temp. Cooling" annotation (Placement(transformation(
          extent={{-9,-10},{9,10}},
          rotation=0,
          origin={-151,-60})));
    Modelica.Blocks.Logical.LessEqualThreshold Le(threshold=T_cool_endday*86400)
      annotation (Placement(transformation(extent={{-200,-110},{-180,-90}})));
    Modelica.Blocks.Logical.And And
      annotation (Placement(transformation(extent={{-160,-90},{-140,-70}})));
    Modelica.Blocks.Interfaces.RealInput H_flo "Valveposition of floor heating"
      annotation (Placement(transformation(
          extent={{-20,-20},{20,20}},
          rotation=0,
          origin={-230,-20})));
    parameter Buildings.HeatTransfer.Data.GlazingSystems.DoubleClearAir13Clear Win(
        haveExteriorShade=true, shade=Buildings.HeatTransfer.Data.Shades.Generic())
      annotation (Placement(transformation(extent={{132,98},{148,114}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin ToC1
      annotation (Placement(transformation(extent={{140,60},{160,80}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin ToC2
      annotation (Placement(transformation(extent={{130,-122},{150,-102}})));
  equation
    connect(ToHeatFlow.y, Roo.qGai_flow) annotation (Line(points={{-159,20},{
            -120,20},{-120,40},{-7.6,40}}, color={0,0,127}));
    connect(ToHeatFlow.u[1], N_Pers)
      annotation (Line(points={{-182,20},{-230,20}}, color={0,0,127}));
    connect(Tsouth.port, WallSouth.port_a)
      annotation (Line(points={{-110,104},{-110,98}},  color={191,0,0}));
    connect(Tceil.port, Ceiling.port_a)
      annotation (Line(points={{-32,104},{-32,98}},  color={191,0,0}));
    connect(Tnorth.port, WallNorth.port_a)
      annotation (Line(points={{-72,104},{-72,98}},  color={191,0,0}));
    connect(TemTab.y, DeMux.u)
      annotation (Line(points={{-159,134},{-152,134}}, color={0,0,127}));
    connect(WallEast.port_a, Teast.port)
      annotation (Line(points={{10,98},{10,104}},  color={191,0,0}));
    connect(Roo.heaPorAir, TrooSen.port) annotation (Line(points={{13,32},{50,
            32},{50,70},{100,70}},
                               color={191,0,0}));
    connect(WallSouth.port_b, Roo.surf_surBou[1]) annotation (Line(points={{-110,78},
            {-110,74},{10.2,74},{10.2,17.2}}, color={191,0,0}));
    connect(Ceiling.port_b, Roo.surf_surBou[4]) annotation (Line(points={{-32,78},
            {-32,74},{10,74},{10,36},{10.2,36},{10.2,18.4}}, color={191,0,0}));
    connect(TsoilSig.y, Tsoil.T)
      annotation (Line(points={{-31,-90},{-15.6,-90}}, color={0,0,127}));
    connect(Val.port_b, FloHeat.port_a)
      annotation (Line(points={{-10,-60},{-4,-60}}, color={0,127,255}));
    connect(Sou.ports[1], Val.port_a)
      annotation (Line(points={{-38,-60},{-30,-60}}, color={0,127,255}));
    connect(WeaData.weaBus, Wb) annotation (Line(
        points={{70,110},{70,50}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%second",
        index=1,
        extent={{-3,-6},{-3,-6}},
        horizontalAlignment=TextAlignment.Right));
    connect(Roo.uSha[1], H_win) annotation (Line(points={{-7.6,50},{-120,50},{
            -120,60},{-230,60}},
                            color={0,0,127}));
    connect(Tsoil.port, FloHeat.surf_b)
      annotation (Line(points={{2,-90},{10,-90},{10,-70}},   color={191,0,0}));
    connect(Mflo.port_b, SinFlo.ports[1])
      annotation (Line(points={{100,-60},{120,-60}}, color={0,127,255}));
    connect(Prod.u2, Mflo.m_flow)
      annotation (Line(points={{98,-6},{90,-6},{90,-49}}, color={0,0,127}));
    connect(Prod.u1, DeltaTflo.y)
      annotation (Line(points={{98,6},{40,6},{40,1}}, color={0,0,127}));
    connect(Prod.y, Cw.u)
      annotation (Line(points={{121,0},{138,0}}, color={0,0,127}));
    connect(Mflo.port_a, Tret.port_b)
      annotation (Line(points={{80,-60},{56,-60}}, color={0,127,255}));
    connect(FloHeat.port_b, Tret.port_a)
      annotation (Line(points={{16,-60},{36,-60}}, color={0,127,255}));
    connect(Tret.T, DeltaTflo.u2)
      annotation (Line(points={{46,-49},{46,-22}}, color={0,0,127}));
    connect(WallNorth.port_b, Roo.surf_surBou[3]) annotation (Line(points={{-72,78},
            {-72,74},{10.2,74},{10.2,18}}, color={191,0,0}));
    connect(WallEast.port_b, Roo.surf_surBou[2])
      annotation (Line(points={{10,78},{10,17.6},{10.2,17.6}}, color={191,0,0}));
    connect(Roo.surf_surBou[5], FloHeat.surf_a) annotation (Line(points={{10.2,18.8},
            {10,18.8},{10,-50}}, color={191,0,0}));
    connect(Sw.y, Sou.T_in) annotation (Line(points={{-97,-68},{-68,-68},{-68,
            -56},{-60,-56}},
                        color={0,0,127}));
    connect(HcvFlo.TSup, Sw.u3) annotation (Line(points={{-126,-87},{-126,-76},
            {-120,-76}},
                   color={0,0,127}));
    connect(Tcool.y, Sw.u1) annotation (Line(points={{-141.1,-60},{-120,-60}},
                              color={0,0,127}));
    connect(P_flo, Cw.y)
      annotation (Line(points={{186,0},{161,0}}, color={0,0,127}));
    connect(Le.u, Wb.solTim) annotation (Line(points={{-202,-100},{-220,-100},{
            -220,-120},{70,-120},{70,50}},
                                      color={0,0,127}), Text(
        string="%second",
        index=1,
        extent={{-6,3},{-6,3}},
        horizontalAlignment=TextAlignment.Right));
    connect(Le.u, Ge.u) annotation (Line(points={{-202,-100},{-220,-100},{-220,
            -60},{-202,-60}},
                         color={0,0,127}));
    connect(Ge.y, And.u1) annotation (Line(points={{-179,-60},{-170,-60},{-170,
            -80},{-162,-80}},
                         color={255,0,255}));
    connect(Le.y, And.u2)
      annotation (Line(points={{-179,-100},{-170,-100},{-170,-88},{-162,-88}},
                                                       color={255,0,255}));
    connect(And.y, Sw.u2)
      annotation (Line(points={{-139,-80},{-130,-80},{-130,-68},{-120,-68}},
                                                       color={255,0,255}));
    connect(Val.y, H_flo) annotation (Line(points={{-20,-48},{-20,-20},{-230,
            -20}},
          color={0,0,127}));
    connect(DeMux.y5[1], Tsouth.T) annotation (Line(points={{-129,126},{-110,
            126},{-110,121.6}},
                         color={0,0,127}));
    connect(DeMux.y4[1], Tnorth.T) annotation (Line(points={{-129,130},{-72,130},
            {-72,121.6}},
                       color={0,0,127}));
    connect(DeMux.y3[1], Tceil.T) annotation (Line(points={{-129,134},{-32,134},
            {-32,121.6}},
                   color={0,0,127}));
    connect(DeMux.y2[1], Teast.T)
      annotation (Line(points={{-129,138},{10,138},{10,121.6}},
                                                              color={0,0,127}));
    connect(TrooSen.T, ToC1.Kelvin)
      annotation (Line(points={{120,70},{138,70}}, color={0,0,127}));
    connect(ToC2.Kelvin, Sou.T_in) annotation (Line(points={{128,-112},{-68,
            -112},{-68,-56},{-60,-56}},
                                  color={0,0,127}));
    connect(DeltaTflo.u1, Sou.T_in) annotation (Line(points={{34,-22},{34,-40},
            {-68,-40},{-68,-56},{-60,-56}},
                                       color={0,0,127}));
    connect(Roo.weaBus, Wb) annotation (Line(
        points={{31.9,49.9},{45.95,49.9},{45.95,50},{70,50}},
        color={255,204,51},
        thickness=0.5), Text(
        string="%second",
        index=1,
        extent={{6,3},{6,3}},
        horizontalAlignment=TextAlignment.Left));
    connect(ToC1.Celsius, T_room)
      annotation (Line(points={{161,70},{186,70}}, color={0,0,127}));
    connect(ToC2.Celsius, T_sup_flo)
      annotation (Line(points={{151,-112},{186,-112}}, color={0,0,127}));
    connect(HcvFlo.TOut, Wb.TDryBul) annotation (Line(points={{-126,-110},{-126,
            -120},{70,-120},{70,50}}, color={0,0,127}), Text(
        string="%second",
        index=1,
        extent={{-3,-6},{-3,-6}},
        horizontalAlignment=TextAlignment.Right));
    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-240,-140},{200,
              160}})),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-240,-140},{
              200,160}}),
                      graphics={
          Rectangle(
            extent={{-190,150},{30,68}},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Text(
            extent={{-44,146},{58,142}},
            lineColor={0,0,0},
            fillColor={255,170,170},
            fillPattern=FillPattern.None,
            fontSize=14,
            fontName="Arial",
            horizontalAlignment=TextAlignment.Left,
            textStyle={TextStyle.Bold},
            textString="Innenwände und Decke"),
          Rectangle(
            extent={{-224,-36},{-90,-130}},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Text(
            extent={{-176,-42},{-64,-46}},
            lineColor={0,0,0},
            fillColor={255,170,170},
            fillPattern=FillPattern.None,
            fontSize=14,
            fontName="Arial",
            horizontalAlignment=TextAlignment.Left,
            textStyle={TextStyle.Bold},
            textString="Vorlauftemperaturregelung"),
          Rectangle(
            extent={{-80,-36},{160,-130}},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Text(
            extent={{-74,-124},{56,-128}},
            lineColor={0,0,0},
            fillColor={255,170,170},
            fillPattern=FillPattern.None,
            fontSize=14,
            fontName="Arial",
            horizontalAlignment=TextAlignment.Left,
            textStyle={TextStyle.Bold},
            textString="Fußboden und Fußbodenheizung"),
          Rectangle(
            extent={{-80,60},{40,12}},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Text(
            extent={{-76,38},{40,6}},
            lineColor={0,0,0},
            fillColor={255,170,170},
            fillPattern=FillPattern.None,
            fontSize=14,
            fontName="Arial",
            horizontalAlignment=TextAlignment.Left,
            textStyle={TextStyle.Bold},
            textString="Raum mit Außen-
wand und Fenster"),
          Rectangle(
            extent={{44,150},{160,90}},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Text(
            extent={{52,146},{168,142}},
            lineColor={0,0,0},
            fillColor={255,170,170},
            fillPattern=FillPattern.None,
            fontSize=14,
            fontName="Arial",
            horizontalAlignment=TextAlignment.Left,
            textStyle={TextStyle.Bold},
            textString="Wetter und Materialparameter")}),
      Documentation(info="<html>
      <p>
      This example represents a simplified room model that can be heated and cooled by floor heating.
      The default size of the room is 6 x 10 x 3 m and the size of the window facing west is 2 x 6 m. The 
      room has one outer wall (with the) window that is aligned west (11.5cm brick, 20cm insulation and 
      2,5cm gypsum board). The inner walls are made of 2,5cm gypsum, 12cm brick and 2,5cm gypsum board.
      The ceiling and floor are massive constructions. The floor boundary is the the soil temp. that is modeled
      as sine with mean of 9 derees, amplitude of 6 K and frequency of 1 year (Minimum at 1st of march).
      All other (four) boundary temperatures are read from a table, which temp. values were recorded in a actual
      office building.   
      <p>
      The input H_flo and H_win represent control signals of the floor heating supply valves and window blinds, 
      resp.. The supply temp. of the floor heating is determined by heating curve in winter and is constant during 
      summer cooling. Further, using input N_Pers the number of occupants in the room can be  
      adjusted (75+25+0W radiant, convective and latent heat per person). Weather Data is used from Chemnitz. 
      The outputs of the model are the room air (T_air), and floor heating supply 
      temperatures (T_sup_flo) and the corresponding heating powers (P_flo).  The room has a window 
      facing  west, which blinds can be controlled using H_win (0 opened 1 closed).
      </html>"),
      experiment(
        StopTime=2592000,
        __Dymola_fixedstepsize=15,
        __Dymola_Algorithm="Dassl"));
  end RoomFlo;
  annotation (uses(Buildings(version="7.0.0"), Modelica(version="3.2.3")));
end SingleRoom;
