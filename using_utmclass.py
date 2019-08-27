import utm_c

gc=utm_c.utm_('WGS84')
print(gc.utm2ll(704462.4938, 9282139.5482, 48,'S'))
print(gc.ll2utm(-6.4910555327, 106.8489137899))
print("x = ", gc.ll2utm(-6.4910555327, 106.8489137899)[0])
print("y = ", gc.ll2utm(-6.4910555327, 106.8489137899)[1])
print("zone = ", gc.ll2utm(-6.4910555327, 106.8489137899)[2])
print("hemis = ", gc.ll2utm(-6.4910555327, 106.8489137899)[3])
