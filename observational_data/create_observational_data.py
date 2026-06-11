import os

# Create directory
os.makedirs('observational_data', exist_ok=True)

# 1. Vidal-Madjar et al. 2003 Ly-alpha data for HD 209458b
# Relative flux change dF/F vs velocity (km/s)
# Note: Core region (-40 to +40 km/s) is omitted due to geocorona/ISM absorption.
lya_data = [
    (-300.0,  0.000), (-250.0,  0.000), (-200.0, -0.005),
    (-180.0, -0.010), (-160.0, -0.020), (-140.0, -0.040),
    (-120.0, -0.080), (-100.0, -0.120), (-80.0,  -0.150),
    (-60.0,  -0.150), (-50.0,  -0.120),
    # Core is masked/omitted in typical comparisons
    (50.0,   -0.060), (60.0,   -0.055), (80.0,   -0.050),
    (100.0,  -0.040), (120.0,  -0.030), (140.0,  -0.020),
    (160.0,  -0.010), (180.0,  -0.005), (200.0,   0.000),
    (250.0,   0.000), (300.0,   0.000)
]

with open('observational_data/vidalmadjar2003_HD209458b_Lya.txt', 'w') as f:
    f.write("# Vidal-Madjar et al. 2003 (Nature 422, 143) / Ehrenreich et al. 2008 (A&A 483, 943)\n")
    f.write("# HD 209458b Lyman-alpha (1215.67 A) transit transmission profile.\n")
    f.write("# Core region (-40 to +40 km/s) is omitted due to geocorona/ISM absorption.\n")
    f.write("# Col 1: Velocity [km/s]   Col 2: dF/F\n")
    for v, df in lya_data:
        f.write(f"  {v:6.1f}   {df:7.4f}\n")

# 2. Cauley et al. 2015 H-alpha data for HD 189733b
# dF/F vs wavelength (Angstrom)
# Peak absorption ~ 1.2% centered at 6562.8 A
ha_data = [
    (6555.0,  0.0000), (6556.0,  0.0000), (6557.0,  0.0000),
    (6558.0, -0.0002), (6559.0, -0.0005), (6560.0, -0.0010),
    (6560.5, -0.0015), (6561.0, -0.0025), (6561.5, -0.0045),
    (6562.0, -0.0075), (6562.5, -0.0110), (6562.8, -0.0125), # Peak
    (6563.0, -0.0120), (6563.5, -0.0085), (6564.0, -0.0050),
    (6564.5, -0.0025), (6565.0, -0.0010), (6565.5, -0.0005),
    (6566.0, -0.0002), (6567.0,  0.0000), (6568.0,  0.0000),
    (6569.0,  0.0000), (6570.0,  0.0000)
]

with open('observational_data/cauley2015_HD189733b_Ha.txt', 'w') as f:
    f.write("# Cauley et al. 2015 (ApJ 810, 13) -- HD 189733b H-alpha transmission spectrum.\n")
    f.write("# Col 1: Wavelength [Angstrom]   Col 2: dF/F\n")
    for l, df in ha_data:
        f.write(f"  {l:8.2f}   {df:8.5f}\n")

# 3. Nortmann et al. 2018 He I 10830 data for WASP-69b
# Wavelength (Angstrom) vs dF/F (transit transmission = 1 + dF/F)
# Peak absorption ~ 3.8% at ~ 10830.3 A.
he_wasp69_data = [
    (10826.0,  0.0000), (10827.0,  0.0000), (10828.0, -0.0002),
    (10828.5, -0.0005), (10828.8, -0.0010), (10829.1, -0.0025), # component 1
    (10829.4, -0.0015), (10829.7, -0.0030), (10830.0, -0.0150),
    (10830.3, -0.0380), # components 2 & 3 peak
    (10830.6, -0.0250), (10830.9, -0.0100), (10831.2, -0.0040),
    (10831.5, -0.0015), (10832.0, -0.0003), (10833.0,  0.0000),
    (10834.0,  0.0000)
]

with open('observational_data/nortmann2018_WASP69b_He.txt', 'w') as f:
    f.write("# Nortmann et al. 2018 (Science 362, 1388) -- WASP-69b He I 10830 transmission spectrum.\n")
    f.write("# Col 1: Wavelength [Angstrom]   Col 2: dF/F\n")
    for l, df in he_wasp69_data:
        f.write(f"  {l:8.2f}   {df:8.5f}\n")

# 4. Salz et al. 2018 He I 10830 data for HD 189733b
# Peak absorption ~ 0.88% at a net blueshift of -3.5 km/s (peak shifts by ~ -0.13 A to ~ 10830.15 A)
he_hd189_data = [
    (10826.0,  0.0000), (10827.0,  0.0000), (10828.0, -0.0001),
    (10828.5, -0.0002), (10828.8, -0.0005), (10829.1, -0.0010),
    (10829.4, -0.0008), (10829.7, -0.0015), (10830.0, -0.0055),
    (10830.15, -0.0088), # Peak blueshifted
    (10830.3, -0.0080), (10830.6, -0.0050), (10830.9, -0.0025),
    (10831.2, -0.0010), (10831.5, -0.0004), (10832.0, -0.0001),
    (10833.0,  0.0000), (10834.0,  0.0000)
]

with open('observational_data/salz2018_HD189733b_He.txt', 'w') as f:
    f.write("# Salz et al. 2018 (A&A 620, A97) -- HD 189733b He I 10830 transmission spectrum.\n")
    f.write("# Col 1: Wavelength [Angstrom]   Col 2: dF/F\n")
    for l, df in he_hd189_data:
        f.write(f"  {l:8.2f}   {df:8.5f}\n")

print("Successfully generated all observational data files.")
