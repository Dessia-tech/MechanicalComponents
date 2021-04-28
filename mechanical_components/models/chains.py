
import mechanical_components.chains as mcc

# https://www.tridistribution.fr/chaines-a-rouleaux-norme-europeene-iso/6184-5453-rollerchain-simplex-european-series.html#/4017-ref-04b1/4018-p-6/17601-w-280/17602-o_d-400/17603-o_d-185/17604-c-830
iso_chain1 = mcc.RollerChain(pitch=0.006, inner_width=0.002,
                          overall_width=0.0083, roller_outer_diameter=0.004,
                          pin_diameter=0.00185)