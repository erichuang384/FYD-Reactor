using Clapeyron
model = CompositeModel(["water","paracetamol"];fluid=PCSAFT,solid=SolidHfus)
T = 298.15
p = 1e5
z = [1.,0.]
s = sle_solubility(model,p,T,[1.,1.];solute=["paracetamol"])

isobaric_heat_capacity(model,p,T,z)