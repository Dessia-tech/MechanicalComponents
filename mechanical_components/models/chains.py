
import mechanical_components.chains as mcc

# https://www.tridistribution.fr/chaines-a-rouleaux-norme-europeene-iso/6184-5453-rollerchain-simplex-european-series.html#/4017-ref-04b1/4018-p-6/17601-w-280/17602-o_d-400/17603-o_d-185/17604-c-830

# https://www.bege.nl/downloads/catalogues/BEGE%20Drive%20Components%202014%20EN.pdf  page 20
iso_chain1 = mcc.RollerChain(pitch=0.006, inner_plates_width=0.0028,
                             outer_plates_width=0.0041,
                             overall_width=0.0066, roller_diameter=0.004,
                             plate_height=0.005,
                             pin_diameter=0.00185)
iso_chain2 = mcc.RollerChain(pitch=0.008, inner_plates_width=0.003,
                             outer_plates_width=0.00477,
                             overall_width=0.0078, roller_diameter=0.005,
                             plate_height=0.00675,
                             pin_diameter=0.00231)

iso_chain3 = mcc.RollerChain(pitch=0.00952, inner_plates_width=0.00572,
                             outer_plates_width=0.00853,
                             overall_width=0.0130, roller_diameter=0.00635,
                             plate_height=0.00826,
                             pin_diameter=0.00328)

iso_chain8 = mcc.RollerChain(pitch=0.01587, inner_plates_width=0.00640,
                             outer_plates_width=0.0108,
                             overall_width=0.0162, roller_diameter=0.01016,
                             plate_height=0.01470,
                             pin_diameter=0.00508)

iso_chains = [iso_chain1, iso_chain2, iso_chain3, iso_chain8]