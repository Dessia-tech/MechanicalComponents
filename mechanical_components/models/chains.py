
import mechanical_components.chains as mcc

# https://www.tridistribution.fr/chaines-a-rouleaux-norme-europeene-iso/6184-5453-rollerchain-simplex-european-series.html#/4017-ref-04b1/4018-p-6/17601-w-280/17602-o_d-400/17603-o_d-185/17604-c-830

# https://www.bege.nl/downloads/catalogues/BEGE%20Drive%20Components%202014%20EN.pdf  page 20
# https://www.fbchain.com/fbkc-roller-chain-solutions-quick-find
# iso_chain1 = mcc.RollerChain(pitch=0.006, inner_plates_width=0.0028,
#                              outer_plates_width=0.0041,
#                              overall_width=0.0066, roller_diameter=0.004,
#                              plate_height=0.005,
#                              pin_diameter=0.00185)
# iso_chain2 = mcc.RollerChain(pitch=0.008, inner_plates_width=0.003,
#                              outer_plates_width=0.00477,
#                              overall_width=0.0078, roller_diameter=0.005,
#                              plate_height=0.00675,
#                              pin_diameter=0.00231)

# iso_chain3 = mcc.RollerChain(pitch=0.00952, inner_plates_width=0.00572,
#                              outer_plates_width=0.00853,
#                              overall_width=0.0130, roller_diameter=0.00635,
#                              plate_height=0.00826,
#                              pin_diameter=0.00328)

# iso_chain8 = mcc.RollerChain(pitch=0.01587, inner_plates_width=0.00640,
#                              outer_plates_width=0.0108,
#                              overall_width=0.0162, roller_diameter=0.01016,
#                              plate_height=0.01470,
#                              pin_diameter=0.00508)

# iso_chains = [iso_chain1, iso_chain2, iso_chain3, iso_chain8]
# http://www.farnell.com/datasheets/17402.pdf
iso_chain_06B1 = mcc.RollerChain(pitch=0.00952,
                                 roller_width=0.00572,
                                 inner_plate_width=0.0013,
                                 outer_plate_width=0.0013,
                                 pin_length=0.01410,
                                 roller_diameter=0.00635,
                                 inner_plate_height=0.0082,
                                 outer_plate_height=0.0082,
                                 pin_diameter=0.00328,
                                 name='iso R606 06B1')

iso_chain_08B1 = mcc.RollerChain(pitch=0.0127,
                                 roller_width=0.00851,
                                 inner_plate_width=0.0016,
                                 outer_plate_width=0.0016,
                                 pin_length=0.0182,
                                 roller_diameter=0.00851,
                                 inner_plate_height=0.01180,
                                 outer_plate_height=0.01180,
                                 pin_diameter=0.00445,
                                 name='iso R606 08B1')

iso_chain_10B1 = mcc.RollerChain(pitch=0.015875,
                                 roller_width=0.00965,
                                 inner_plate_width=0.0017,
                                 outer_plate_width=0.0017,
                                 pin_length=0.02090,
                                 roller_diameter=0.01016,
                                 inner_plate_height=0.0147,
                                 outer_plate_height=0.0147,
                                 pin_diameter=0.00508,
                                 name='iso R606 10B1')


iso_chain_12B1 = mcc.RollerChain(pitch=0.01905,
                                 roller_width=0.01168,
                                 inner_plate_width=0.00185,
                                 outer_plate_width=0.00185,
                                 pin_length=0.02420,
                                 roller_diameter=0.01207,
                                 inner_plate_height=0.016,
                                 outer_plate_height=0.016,
                                 pin_diameter=0.00572,
                                 name='iso R606 12B1')

iso_chain_16B1 = mcc.RollerChain(pitch=0.0254,
                                 roller_width=0.01702,
                                 inner_plate_width=0.00415,
                                 outer_plate_width=0.0031,
                                 pin_length=0.0374,
                                 roller_diameter=0.01588,
                                 inner_plate_height=0.021,
                                 outer_plate_height=0.021,
                                 pin_diameter=0.00828,
                                 name='iso R606 16B1')

iso_chains = [iso_chain_06B1, iso_chain_08B1, iso_chain_10B1, iso_chain_12B1,
              iso_chain_16B1]