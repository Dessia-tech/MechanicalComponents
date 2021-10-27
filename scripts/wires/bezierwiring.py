import volmdlr as vm
import volmdlr.edges as edges
import volmdlr.primitives3d as p3d
import mechanical_components.wires as wires
import math
import volmdlr.wires as vmw

route_ways = [[vm.Point3D(-0.5089699401855469, -0.04520467758178711,
                          0.36238052368164064),
               vm.Point3D(-0.5053260803222657, -0.04520467758178711,
                          0.35069119262695314),
               vm.Point3D(-0.4948938903808594, -0.04520467758178711,
                          0.34381967163085936),
               vm.Point3D(-0.48224090576171874, -0.07058212280273438,
                          0.3488803405761719),
               vm.Point3D(-0.45849005126953124, -0.11041954040527344,
                          0.3169632263183594),
               vm.Point3D(-0.420466796875, -0.15144271850585939,
                          0.2875592041015625),
               vm.Point3D(-0.44445648193359377, -0.1814255828857422,
                          0.24849771118164063)],
              [vm.Point3D(-0.5018493957519531, 0.18525401306152345,
                          0.35413470458984375),
               vm.Point3D(-0.49670852661132814, 0.23336582946777346,
                          0.32308316040039065),
               vm.Point3D(-0.49799542236328126, 0.24286868286132812,
                          0.33193890380859375),
               vm.Point3D(-0.498003662109375, 0.24293011474609374,
                          0.33199615478515626),
               vm.Point3D(-0.4986155700683594, 0.25198350524902347,
                          0.3365279541015625),
               vm.Point3D(-0.49862768554687503, 0.25309474182128905,
                          0.336674072265625),
               vm.Point3D(-0.4906199645996094, 0.2645937194824219,
                          0.34445428466796874),
               vm.Point3D(-0.48758999633789063, 0.26885696411132814,
                          0.3473858947753906),
               vm.Point3D(-0.47301531982421874, 0.27804080200195314,
                          0.359903564453125),
               vm.Point3D(-0.46485189819335937, 0.2814581604003906,
                          0.36667327880859374),
               vm.Point3D(-0.46122686767578125, 0.2800146179199219,
                          0.36667327880859374),
               vm.Point3D(-0.45583, 0.28778, 0.36882)],
              [vm.Point3D(-0.5018493957519531, 0.18525401306152345,
                          0.35413470458984375),
               vm.Point3D(-0.5047180480957031, 0.24286868286132812,
                          0.37570556640625),
               vm.Point3D(-0.5094658813476562, 0.28496697998046877,
                          0.4091636047363281),
               vm.Point3D(-0.5094107666015625, 0.29067028808593753,
                          0.4091499328613281),
               vm.Point3D(-0.5086275329589844, 0.29269387817382814,
                          0.4115304870605469),
               vm.Point3D(-0.492208984375, 0.34778594970703125,
                          0.41672152709960936),
               vm.Point3D(-0.49219058227539064, 0.34889974975585936,
                          0.4167613525390625),
               vm.Point3D(-0.4906031188964844, 0.3503307800292969,
                          0.416326171875),
               vm.Point3D(-0.43997076416015624, 0.3820563049316406,
                          0.41754165649414066)],
              [vm.Point3D(-0.48461538696289064, 0.07002466583251954,
                          0.3288733825683594),
               vm.Point3D(-0.47543768310546874, 0.1276393356323242,
                          0.2743473815917969),
               vm.Point3D(-0.45788131713867186, 0.18525401306152345,
                          0.2596198425292969),
               vm.Point3D(-0.4491347351074219, 0.18525401306152345,
                          0.25364218139648437)],
              [vm.Point3D(-0.44445648193359377, -0.1814255828857422,
                          0.24849771118164063),
               vm.Point3D(-0.44445648193359377, -0.220064697265625,
                          0.406675537109375),
               vm.Point3D(-0.38691531372070315, -0.21795477294921875,
                          0.4165601501464844),
               vm.Point3D(-0.3293741455078125, -0.25439987182617185,
                          0.42312609863281253),
               vm.Point3D(-0.2718329772949219, -0.27062741088867187,
                          0.42724597167968753),
               vm.Point3D(-0.21429180908203124, -0.2835421142578125,
                          0.43037698364257815)],
              [vm.Point3D(-0.43997076416015624, 0.3820563049316406,
                          0.41754165649414066),
               vm.Point3D(-0.3813824768066406, 0.4160341491699219,
                          0.4684100036621094),
               vm.Point3D(-0.3420126037597656, 0.39170513916015626,
                          0.4684100036621094),
               vm.Point3D(-0.28899182128906253, 0.3953202514648438,
                          0.49710830688476565),
               vm.Point3D(-0.2232827911376953, 0.393857177734375,
                          0.49710830688476565)],
              [vm.Point3D(-0.43997076416015624, 0.3820563049316406,
                          0.41754165649414066),
               vm.Point3D(-0.4444637451171875, 0.42065655517578127,
                          0.36667327880859374),
               vm.Point3D(-0.41602664184570315, 0.3660671997070313,
                          0.36667327880859374)],
              [vm.Point3D(-0.4491347351074219, 0.18525401306152345,
                          0.25364218139648437),
               vm.Point3D(-0.43980322265625, 0.18525401306152345,
                          0.24938922119140625),
               vm.Point3D(-0.4296233520507813, 0.18525401306152345,
                          0.2468501434326172),
               vm.Point3D(-0.4111953735351563, 0.18525401306152345,
                          0.24407789611816408),
               vm.Point3D(-0.4040201110839844, 0.183547607421875,
                          0.2433340301513672),
               vm.Point3D(-0.34729000000000004, 0.16754, 0.30205000000000004)],
              [vm.Point3D(-0.4491347351074219, 0.18525401306152345,
                          0.25364218139648437),
               vm.Point3D(-0.44449313354492187, 0.20157859802246095,
                          0.25558018493652346),
               vm.Point3D(-0.4415464782714844, 0.20164353942871094,
                          0.2140682067871094),
               vm.Point3D(-0.4195151977539063, 0.20167825317382812,
                          0.17792898559570314),
               vm.Point3D(-0.4259250793457031, 0.20165774536132813,
                          0.14754739379882814),
               vm.Point3D(-0.43073440551757813, 0.201645263671875,
                          0.1365624237060547),
               vm.Point3D(-0.4431706237792969, 0.20164329528808594,
                          0.11233148193359375),
               vm.Point3D(-0.44005331420898436, 0.2090709228515625,
                          0.11233148193359375)],
              [vm.Point3D(-0.44445648193359377, -0.17054653930664063,
                          0.1792205505371094),
               vm.Point3D(-0.38691531372070315, -0.15818742370605468,
                          0.07942154693603516),
               vm.Point3D(-0.38691531372070315, -0.19686361694335938,
                          0.042002117156982426),
               vm.Point3D(-0.38691531372070315, -0.20636404418945312,
                          0.03787649536132812),
               vm.Point3D(-0.38691531372070315, -0.21750837707519532,
                          0.037910995483398435),
               vm.Point3D(-0.38691531372070315, -0.22649922180175783,
                          0.04229977798461914),
               vm.Point3D(-0.38691531372070315, -0.264634521484375,
                          0.07554747009277343)],
              [vm.Point3D(-0.2232827911376953, 0.393857177734375,
                          0.49710830688476565),
               vm.Point3D(-0.2170585174560547, 0.39243637084960936,
                          0.5670205078125),
               vm.Point3D(-0.270635986328125, 0.34786636352539063,
                          0.6369326782226563),
               vm.Point3D(-0.24387889099121093, 0.33339053344726566,
                          0.6369326782226563),
               vm.Point3D(-0.23666970825195313, 0.3226582946777344,
                          0.6369326782226563),
               vm.Point3D(-0.23632611083984376, 0.3220267944335938,
                          0.6369326782226563),
               vm.Point3D(-0.22591857910156252, 0.30289871215820313,
                          0.6369326782226563),
               vm.Point3D(-0.22469691467285158, 0.30065341186523437,
                          0.6369326782226563)],
              [vm.Point3D(-0.41602664184570315, 0.3660671997070313,
                          0.36667327880859374),
               vm.Point3D(-0.3714701538085938, 0.35304376220703126,
                          0.33071743774414064),
               vm.Point3D(-0.3556944580078125, 0.36206179809570316,
                          0.315804931640625),
               vm.Point3D(-0.3325960998535156, 0.3750653381347656,
                          0.28737176513671875),
               vm.Point3D(-0.34811000000000003, 0.40676, 0.29207)],
              [vm.Point3D(-0.41602664184570315, 0.3660671997070313,
                          0.36667327880859374),
               vm.Point3D(-0.3819021301269531, 0.3005599670410156,
                          0.36667327880859374),
               vm.Point3D(-0.3721069946289063, 0.27992718505859376,
                          0.38934103393554687),
               vm.Point3D(-0.3686797485351563, 0.2779041748046875,
                          0.3923334655761719),
               vm.Point3D(-0.38975, 0.25953, 0.42764)],
              [vm.Point3D(-0.38691531372070315, -0.264634521484375,
                          0.07554747009277343),
               vm.Point3D(-0.3293741455078125, -0.3256164855957031,
                          0.04978006362915039),
               vm.Point3D(-0.33951220703125, -0.3206627807617188,
                          0.0702730941772461),
               vm.Point3D(-0.36225, -0.31244, 0.05258)],
              [vm.Point3D(-0.38691531372070315, -0.264634521484375,
                          0.07554747009277343),
               vm.Point3D(-0.38691531372070315, -0.3027698059082031,
                          0.108795166015625),
               vm.Point3D(-0.3786444091796875, -0.30464111328125,
                          0.11551530456542969),
               vm.Point3D(-0.38325143432617187, -0.31129998779296875,
                          0.1122625961303711),
               vm.Point3D(-0.38325143432617187, -0.32617294311523437,
                          0.11109738159179687),
               vm.Point3D(-0.38325143432617187, -0.35732135009765625,
                          0.11190681457519532)],
              [vm.Point3D(-0.15675062561035158, -0.2600059814453125,
                          0.38612850952148436),
               vm.Point3D(-0.09920946502685547, -0.2010964813232422,
                          0.34828054809570314),
               vm.Point3D(-0.12727, -0.19047, 0.34784)],
              [vm.Point3D(-0.14443196105957032, 0.39210150146484374,
                          0.49710830688476565),
               vm.Point3D(-0.07504000091552734, 0.39370001220703127,
                          0.5670205078125),
               vm.Point3D(-0.0637728500366211, 0.38746368408203125,
                          0.5670205078125),
               vm.Point3D(-0.06164128494262695, 0.37732806396484375,
                          0.5670205078125),
               vm.Point3D(-0.06123465347290039, 0.3661609497070313,
                          0.5670205078125),
               vm.Point3D(-0.062300304412841796, 0.351792236328125,
                          0.5670205078125),
               vm.Point3D(-0.07594982147216797, 0.32272637939453125,
                          0.5804618530273438)],
              [vm.Point3D(-0.11940857696533204, 0.09856647491455078,
                          0.6369326782226563),
               vm.Point3D(-0.10414730072021484, 0.06612438201904297,
                          0.614640625),
               vm.Point3D(-0.10319176483154296, 0.06254248809814453,
                          0.6133880004882812),
               vm.Point3D(-0.10107452392578126, 0.045622203826904295,
                          0.6044317016601563),
               vm.Point3D(-0.09986939239501953, 0.03692203903198242,
                          0.5995250854492188),
               vm.Point3D(-0.09929310607910156, 0.03290122985839844,
                          0.5972074584960938),
               vm.Point3D(-0.09784471130371095, 0.026640857696533203,
                          0.59435693359375),
               vm.Point3D(-0.0968069839477539, 0.02187256622314453,
                          0.5923446655273438),
               vm.Point3D(-0.09731339263916015, 0.019717880249023437,
                          0.5898930053710938),
               vm.Point3D(-0.09809609222412109, 0.01961002349853516,
                          0.5885726318359376),
               vm.Point3D(-0.09876129913330078, 0.019472860336303712,
                          0.5875516967773438),
               vm.Point3D(-0.09821088409423828, 0.0163995418548584,
                          0.585931396484375)],
              [vm.Point3D(-0.45039691162109374, 0.23704273986816407,
                          0.07763520812988281),
               vm.Point3D(-0.4258146667480469, 0.2941393127441406,
                          0.07763520812988281),
               vm.Point3D(-0.4053294372558594, 0.34171978759765625,
                          0.07763520812988281),
               vm.Point3D(-0.38408477783203127, 0.3660409240722656,
                          0.07763520812988281),
               vm.Point3D(-0.32280224609375, 0.370546630859375,
                          0.07763520812988281),
               vm.Point3D(-0.3112119445800781, 0.3712756042480469,
                          0.07763520812988281),
               vm.Point3D(-0.2991233825683594, 0.3716849670410156,
                          0.07763520812988281),
               vm.Point3D(-0.2552846527099609, 0.37280926513671875,
                          0.07763520812988281),
               vm.Point3D(-0.23685884094238283, 0.3725736083984375,
                          0.07763520812988281),
               vm.Point3D(-0.21114625549316407, 0.3696173095703125,
                          0.07763520812988281),
               vm.Point3D(-0.125942626953125, 0.3584390869140625,
                          0.14248915100097656),
               vm.Point3D(-0.125942626953125, 0.35958511352539063,
                          0.14497906494140625),
               vm.Point3D(-0.125942626953125, 0.36413375854492186,
                          0.14754739379882814),
               vm.Point3D(-0.125942626953125, 0.3781520080566406,
                          0.18371435546875),
               vm.Point3D(-0.10498, 0.35442, 0.17753)],
              [vm.Point3D(-0.38325143432617187, -0.35732135009765625,
                          0.11190681457519532),
               vm.Point3D(-0.33616864013671877, -0.40366592407226565,
                          0.11721980285644532),
               vm.Point3D(-0.28908581542968753, -0.4314606323242188,
                          0.14201609802246093),
               vm.Point3D(-0.30597, -0.4359, 0.15737)],
              [vm.Point3D(-0.07594982147216797, 0.32272637939453125,
                          0.5804618530273438),
               vm.Point3D(-0.08256660461425781, 0.31685479736328126,
                          0.5877514038085938),
               vm.Point3D(-0.0939497299194336, 0.2812198486328125,
                          0.5978882446289062)],
              [vm.Point3D(-0.09821088409423828, 0.0163995418548584,
                          0.585931396484375),
               vm.Point3D(-0.09350689697265625, -0.011791036605834962,
                          0.5717717895507812),
               vm.Point3D(-0.0922275390625, -0.01945809555053711,
                          0.5679207763671875),
               vm.Point3D(-0.09192844390869141, -0.021250526428222657,
                          0.5670205078125),
               vm.Point3D(-0.09145029449462891, -0.03948785400390625,
                          0.5538040161132812),
               vm.Point3D(0.020829999999999998, 0.00218, 0.44610000000000005)],
              [vm.Point3D(-0.0939497299194336, 0.2812198486328125,
                          0.5978882446289062),
               vm.Point3D(-0.10665657806396485, 0.23486692810058593,
                          0.614640625),
               vm.Point3D(-0.10585691070556641, 0.22051191711425783,
                          0.614640625),
               vm.Point3D(-0.10705780029296875, 0.20539178466796876,
                          0.614640625),
               vm.Point3D(-0.10759265899658203, 0.19865748596191407,
                          0.614640625),
               vm.Point3D(-0.10650088500976562, 0.17905870056152345,
                          0.614640625)],
              [vm.Point3D(-0.10650088500976562, 0.17905870056152345,
                          0.614640625),
               vm.Point3D(-0.104323486328125, 0.1399715576171875, 0.614640625),
               vm.Point3D(-0.10663780975341797, 0.11083231353759766,
                          0.614640625),
               vm.Point3D(-0.10414730072021484, 0.06612438201904297,
                          0.614640625),
               vm.Point3D(-0.10319176483154296, 0.06254248809814453,
                          0.6133880004882812),
               vm.Point3D(-0.10107452392578126, 0.045622203826904295,
                          0.6044317016601563),
               vm.Point3D(-0.09986939239501953, 0.03692203903198242,
                          0.5995250854492188),
               vm.Point3D(-0.09929310607910156, 0.03290122985839844,
                          0.5972074584960938),
               vm.Point3D(-0.09784471130371095, 0.026640857696533203,
                          0.59435693359375),
               vm.Point3D(-0.0968069839477539, 0.02187256622314453,
                          0.5923446655273438),
               vm.Point3D(-0.09731339263916015, 0.019717880249023437,
                          0.5898930053710938),
               vm.Point3D(-0.09809609222412109, 0.01961002349853516,
                          0.5885726318359376),
               vm.Point3D(-0.09876129913330078, 0.019472860336303712,
                          0.5875516967773438),
               vm.Point3D(-0.09821088409423828, 0.0163995418548584,
                          0.585931396484375)],
              [vm.Point3D(-0.26766458129882814, -0.012839525222778321,
                          0.6369326782226563),
               vm.Point3D(-0.26869094848632813, -0.016203348159790038,
                          0.6293789672851563),
               vm.Point3D(-0.26986126708984376, -0.020038923263549804,
                          0.6207659301757813),
               vm.Point3D(-0.26769189453125003, -0.02336910629272461,
                          0.6176869506835938),
               vm.Point3D(-0.25792000000000004, -0.0265, 0.6078)]]

new_diameters = [0.020779733838059502, 0.017, 0.02637755832643195,
                 0.010071067811865475, 0.015488088481701516,
                 0.02512461179749811, 0.01021110255092798, 0.007,
                 0.0088309518948453,
                 0.015677078252031313, 0.017083045973594575,
                 0.009000000000000001, 0.007, 0.01, 0.01106225774829855, 0.007,
                 0.019317821063276355, 0.012219544457292886, 0.007,
                 0.01, 0.0174899959967968, 0.009000000000000001,
                 0.0153440804327886, 0.012219544457292886, 0.007]
# b_curves = []
# for i, points in enumerate(route_ways[0]):
#     if len(points) < 3:
#         list_points = []
#         for point1, point2 in zip(points[:-1], points[1:]):
#             if point1 not in list_points:
#                 list_points.append(point1)
#             middle_point = (point1 + point2)*0.5
#             list_points.append(middle_point)
#             if point2 not in list_points:
#                 list_points.append(point2) 
#     else:
#         list_points = points
#     if len(list_points) == 3:
#         bc = edges.BezierCurve3D(2, list_points)
#     else:
#         bc = edges.BezierCurve3D(3, list_points)


#     ax=bc.plot(color = 'b')
#     for point1, point2 in zip(points[:-1], points[1:]):
#         edges.LineSegment3D(point1, point2).plot(ax=ax, color = 'r')

#     # print('3x Dimater :', new_diameters[i]*3)
#     # distances = [point1.point_distance(point2) for point1, point2 in zip(points[:-1], points[1:])]
#     # print('distances :', distances)
#     b_curves.append(bc)


# =============================================================================

# list_junc_wires = []
# for i, points in enumerate([route_ways[1]]):
#     # points = points[::-1]
#     list_points = [points[0]]
#     for i, p in enumerate(points[1:-1]):
#         vec1 = p - list_points[-1]
#         vec2 = points[i+2] - p
#         angle = vm.core.vectors3d_angle(vec1, vec2)
#         # print('angle :', angle*180/math.pi)
#         if angle*180/math.pi > 20:
#             list_points.append(p)
#     # print('distances:', [point1.point_distance(point2) for point1, point2 in zip(points[:-1], points[1:])])
#     list_points.append(points[-1])
#     points = list_points
#     # 
#     # print('list_points :', list_points)
#     list_points = [points[0]]
#     for i, p in enumerate(points[1:]):
#         distance = list_points[-1].point_distance(p)
#         if distance > 2e-3:
#             list_points.append(p)
            
#     # list_points = [points[0]]
#     # for p1, p2 in zip(points[1:-1][1:], points[1:-1][:-1]):
#     #     mid_point = (p1+p2)*0.5
#     #     list_points.append(mid_point)
#     # list_points.append(points[-1])
#     points = list_points
#     # print('len route, len list points :', (len(route_ways[1]), len(list_points)))
#     # print('points:', points)
#     # points.remove(points[-3])
#     diameter = new_diameters[1]
#     # j_w = wires.JunctionWire.n_junction_wires(points, diameter*3, diameter)
#     # list_junc_wires.append(j_w)
#     # # print('distances:', [point1.point_distance(point2) for point1, point2 in zip(points[:-1], points[1:])])
#     # if len(points) == 3:
#     #     bc = edges.BezierCurve3D(2, points)
#     # else:
#     #     bc = edges.BezierCurve3D(3, points)
    
#     # ax=bc.plot(color = 'b')
#     # for point1, point2 in zip(points[:-1], points[1:]):
#     #     edges.LineSegment3D(point1, point2).plot(ax=ax, color = 'r')
#     # for point in points:
#     #     point.plot(ax=ax)
#     # j_w.plot(ax=ax)

# # =============================================================================
# # for i, point in enumerate(points[:-1]):
# list_wires = []
# p1,p2, p3 = points[0], points[1], points[-1]
# t1, t3 = points[1] - points[0], points[-2] - points[-1]

# if len(points) < 3:
#     orls3d = p3d.OpenRoundedLineSegments3D(points, {})
#     list_wires = [orls3d]

# elif len(points) == 3:
#     t2 = -t1 + t3
#     t1.normalize()
#     t2.normalize()
#     t3.normalize()
    
# elif len(points) > 3:
#     orls3d = p3d.OpenRoundedLineSegments3D(points, {})
#     length = 0
#     # for prim in orls3d.primitives[1:-1]:
#     #     if prim.length() > length:
#     #         length = prim.length()
#     #         p2 = (prim.start + prim.end)*0.5
#     #         t2 = prim.start - prim.end
            
#     t1.normalize()
#     t3.normalize()
#     index = int(len(points)/2) - 1
#     p2 = (points[index] + points[index+1])*0.5
#     t2 = points[index] - points[index+1]
#     t2.normalize()
    
# new_points = [p1, p2, p3]
# tangents = [t1, t2, t3]
# tangent1 = t1
# for i in range(2):
    
#     lmin = new_points[i].point_distance(new_points[i+1])
#     lmax = 3*lmin
#     # try:
#     wire = wires.JunctionWire.curvature_radius(new_points[i], tangent1, new_points[i+1], tangents[i+1], 3*diameter, diameter, lmin, lmax)
#     #     valid = True
#     # except ValueError:
#     #     valid = False
#     # while not valid:
#     #     lmax += 0.001
#     #     try:
#     #         wire = wires.JunctionWire.curvature_radius(new_points[i], tangent1, new_points[i+1], tangents[i+1], 3*diameter, diameter, lmin, lmax)
#     #         valid = True
#     #     except ValueError:
#     #         valid = False
#     if wire.path.minimum_radius()[0] < 3*diameter:
#         junction_wire = wires.JunctionWireTangentOptimizer(junction_wire=wire, tangent2 = tangents[i+1],
#                                                            range_theta = [0.0,10.0], range_alpha = [0.0,360.0], Lmax = lmax, Lmin=lmin).optimize_tangent()
#         tangent2 =  junction_wire.tangeancy2
#     list_wires.append(wire)
#     tangent1 = -tangents[i+1]

# list_wires = [list_wires[0].path, list_wires[1].path]
# wire = vmw.Wire3D(list_wires)
# if len(points) == 3:
#     bc = edges.BezierCurve3D(2, points)
# else:
#     bc = edges.BezierCurve3D(3, points)

# ax=bc.plot(color = 'b')
# for point1, point2 in zip(points[:-1], points[1:]):
#     edges.LineSegment3D(point1, point2).plot(ax=ax, color = 'r')
# for point in points:
#     point.plot(ax=ax)
# wire.plot(ax=ax)

# =============================================================================
# Test Other method
# =============================================================================


from scipy.optimize import minimize, least_squares

def find_interval_fixation(orls3d, fixation_rule_from_start_end, dist_max_between_fix):
    #à partir de quel abscissa le point est à fixation_rule_from_start_end du début
    length = orls3d.length()
    
    def find_u_min_from_start(x):
        point = orls3d.point_at_abscissa(x[0]*length)
        
        d1 = abs(orls3d.points[0].point_distance(point)-fixation_rule_from_start_end)
        d2 = abs(x[0]*length-fixation_rule_from_start_end)
        return 100*d1 + 100*d2
        # return abs(orls3d.points[0].point_distance(point)-fixation_rule_from_start_end)*1000
    
    def find_u_max_from_start(x):
        point = orls3d.point_at_abscissa(x[0]*length)
        d1 = abs(orls3d.points[0].point_distance(point)-(dist_max_between_fix+fixation_rule_from_start_end))
        d2 = abs(x[0]*length-(dist_max_between_fix+fixation_rule_from_start_end))
        return 100*d1 + 100*d2
        # return abs(orls3d.points[0].point_distance(point)-(dist_max_between_fix+fixation_rule_from_start_end))*1000
    
    def find_u_min_from_end(x):
        point = orls3d.point_at_abscissa(x[0]*length)
        d1 = abs(orls3d.points[-1].point_distance(point)-fixation_rule_from_start_end)
        d2 = abs((1-x[0])*length-fixation_rule_from_start_end)
        return 100*d1 + 100*d2
        # return abs(orls3d.points[-1].point_distance(point)-fixation_rule_from_start_end)*1000
    
    def find_u_max_from_end(x):
        point = orls3d.point_at_abscissa(x[0]*length)
        d1 = abs(orls3d.points[-1].point_distance(point)-(dist_max_between_fix+fixation_rule_from_start_end))
        d2 = abs((1-x[0])*length-(dist_max_between_fix+fixation_rule_from_start_end))
        return 100*d1 + 100*d2
        # return abs(orls3d.points[-1].point_distance(point)-(dist_max_between_fix+fixation_rule_from_start_end))*1000
    
    
    min1 = least_squares(find_u_min_from_start, (0.01), bounds=(0,1))
    max1 = least_squares(find_u_max_from_start, (0.99), bounds=(0,1))
    min2 = least_squares(find_u_min_from_end, (0.99), bounds=(0,1))
    max2 = least_squares(find_u_max_from_end, (0.01), bounds=(0,1))
    
    minmax1_minmax2 = [] 
    for sol in [min1.x[0], max1.x[0], min2.x[0], max2.x[0]]:
        if math.isclose(sol,0, abs_tol=1e-6):
            minmax1_minmax2.append(0)
        elif math.isclose(sol,1, abs_tol=1e-6):
            minmax1_minmax2.append(1)
        else :
            minmax1_minmax2.append(sol)
        
    return minmax1_minmax2 


fixation_rule_from_start_end = 40e-3
dist_max_between_fix = 200e-3

# diam = new_diameters[1]
# targeted_curv = 3*diam
list_junc_wires = []
for diam, points in zip([new_diameters[13]],[route_ways[13]]):#zip([new_diameters[-2]],[route_ways[-2]]):
    targeted_curv = 3*diam
    
    orls3d = p3d.OpenRoundedLineSegments3D(points, {})
    length = orls3d.length()
    print('length',length)
    
    if len(points)<3 :
        jw = [orls3d]
    
    elif length <= 2*fixation_rule_from_start_end :
        
        tan1 = points[1]-points[0]
        tan1.normalize()
        
        tan2 = points[-2]-points[-1]
        tan2.normalize()
        
        lmin = points[0].point_distance(points[-1])
        jw = [wires.JunctionWire.curvature_radius(point1 = points[0], 
                                                  tangeancy1 = tan1, 
                                                  point2 = points[-1],
                                                  tangeancy2 = tan2, 
                                                  targeted_curv = targeted_curv, 
                                                  diameter = diam, 
                                                  length_min = lmin, 
                                                  length_max = 3*lmin)]
            
    elif length > 2*fixation_rule_from_start_end and \
        length <= 2*fixation_rule_from_start_end + dist_max_between_fix :
            
        ####Solution with fixation
        min1, max1, min2, max2 =  find_interval_fixation(orls3d, 
                                                          fixation_rule_from_start_end,
                                                          dist_max_between_fix)
        
        fix_interval = []
        int1, int2 = [min1, max1], [min2, max2]
        
        ax=points[0].plot()
        for pt in points[1:]:
            pt.plot(ax=ax)
        for point1, point2 in zip(points[:-1], points[1:]):
            edges.LineSegment3D(point1, point2).plot(ax=ax, color = 'r')
            
        point1 = orls3d.point_at_abscissa(min1*length)
        point1.plot(ax=ax, color = 'm')
        point2 = orls3d.point_at_abscissa(max1*length)
        point2.plot(ax=ax, color = 'm')
        point3 = orls3d.point_at_abscissa(min2*length)
        point3.plot(ax=ax, color = 'g')
        point4 = orls3d.point_at_abscissa(max2*length)
        point4.plot(ax=ax, color = 'g')
        
        # # print(point1, point2, point3, point4)
        
        print('distance')
        print(min1, max1, min2, max2)
        print(points[0].point_distance(point1), points[0].point_distance(point2))
        print(points[-1].point_distance(point3), points[-1].point_distance(point4))
        print(points[0].point_distance(points[-1]))
        
        
        if min2 < max1 and min2 > min1:
        #         min1[                      ]max1
        #     0------------------------------1
        # max2[             ]min2
            if max2 < min1 :
                fix_interval.append([min1, min2])
            else :
                fix_interval.append([max2, min2])
                
        elif max2 > min1 and max2 < max1 :
        #         min1[                      ]max1
        #     0------------------------------1
        #              max2[             ]min2
            if min2 < max1 :
                fix_interval.append([max2, min2])
            else :
                fix_interval.append([max2, max1])
            
        point1 = orls3d.point_at_abscissa(fix_interval[0][0]*length)
        point2 = orls3d.point_at_abscissa(fix_interval[0][1]*length)
        
        print('distance')
        print(points[0].point_distance(point1), points[-1].point_distance(point2))
        
        ni = 20
        
        for i, interval in enumerate(fix_interval) :
            if i == 0:
                solutions = []
                for n in range(ni):
                    #try ni position in the interval to find the best fixation
                    abscissa_fix = interval[0] + n*(interval[1]-interval[0])/(ni-1)
                    fix = orls3d.point_at_abscissa(abscissa_fix*length)
                    
                    
                    new_points = [points[0], fix, points[-1]]
                    abscissa_coord = [0, abscissa_fix, 1]
                    
                    tangents, points_jw = [], []
                    for k in range(len(new_points)-1):
                        p1, p2 = new_points[k:k+2]
                        abs1, abs2 = abscissa_coord[k:k+2]
                        
                        tan1 = orls3d.point_at_abscissa((abs1+1/100)*length) - p1
                        tan1.normalize()
                        
                        tan2 = orls3d.point_at_abscissa((abs2-1/100)*length) - p2
                        tan2.normalize()
                        
                        tangents.append([tan1, tan2])
                        
                        points_jw.append([p1, p2])
                    
                    jw = []
                    for pts, tans in zip(points_jw, tangents):
                        lmin = pts[0].point_distance(pts[1])
                        lmax= 2*lmin
                        
                        jw_to_add = wires.JunctionWire.curvature_radius(point1 = pts[0],
                                                                        tangeancy1 = tans[0], 
                                                                        point2 = pts[1],
                                                                        tangeancy2 = tans[1], 
                                                                        targeted_curv = targeted_curv, 
                                                                        diameter = diam, 
                                                                        length_min = lmin, 
                                                                        length_max = lmax)
                        try : 
                            jw_to_add.path
                            jw.append(jw_to_add)
                        except ValueError:
                            # print()
                            # print('lmin', lmin)
                            # print('lmax', lmax)
                            # print('jw_to_add', jw_to_add.targeted_length)
                            # print()
                            break
                        
                    solutions.append(jw)
                    
                good_solutions = []    
                for sol in solutions :
                    radius_curves = []
                    if sol :
                        # ax=points[0].plot()
                        # for pt in points[1:]:
                            # pt.plot(ax=ax)
                        # for point1, point2 in zip(points[:-1], points[1:]):
                            # edges.LineSegment3D(point1, point2).plot(ax=ax, color = 'r')
                            
                        
                        for jw in sol :
                            # jw.path.plot(ax=ax)
                            # for pt in jw.waypoints :
                                # pt.plot(ax=ax, color='g')
                            radius_curves.append(jw.path.minimum_radius())
                    # print(radius_curves, targeted_curv)
                    valid = True
                    for r in radius_curves :
                        if r < targeted_curv :
                            valid = False
                            break
                    if valid and len(radius_curves)>1:
                        good_solutions.append(sol)
                
                # print(len(good_solutions), 'solutions')
                
            ####Solution without fixation if length < dist_max_between_fix
            
            if length < dist_max_between_fix:
                
                p1, p2 = points[0], points[-1]
                
                tan1 = points[1]-p1
                tan1.normalize()
                
                tan2 = points[-2]-p2
                tan2.normalize()
                
                # tan1 = orls3d.point_at_abscissa((1/100)*length) - p1
                # tan1.normalize()
                # tan2 = orls3d.point_at_abscissa((99/100)*length) - p2
                # tan2.normalize()
                
                lmin = p1.point_distance(p2)
                lmax= 2*lmin
                
                jw_to_add = wires.JunctionWire.curvature_radius(point1 = p1, 
                                                                tangeancy1 = tan1, 
                                                                point2 = p2,
                                                                tangeancy2 = tan2, 
                                                                targeted_curv = targeted_curv, 
                                                                diameter = diam, 
                                                                length_min = lmin, 
                                                                length_max = lmax)
                try : 
                    jw_to_add.path
                    
                except ValueError:
                    # print()
                    # print('lmin', lmin)
                    # print('lmax', lmax)
                    # print('jw_to_add', jw_to_add.targeted_length)
                    # print()
                    break
                
                # ax=points[0].plot()
                # for pt in points[1:]:
                #     pt.plot(ax=ax)
                # for point1, point2 in zip(points[:-1], points[1:]):
                #     edges.LineSegment3D(point1, point2).plot(ax=ax, color = 'r')
                        
                    
                # jw_to_add.path.plot(ax=ax)
                # for pt in jw_to_add.waypoints :
                    # pt.plot(ax=ax, color='g')
                radius_curves = jw_to_add.path.minimum_radius()
                # print(radius_curves, targeted_curv)
                if radius_curves > targeted_curv :
                    good_solutions.append(jw_to_add)
                    
            print(len(good_solutions), 'solutions')
            list_junc_wires.append(good_solutions)
            
                
                        
                        
                        
                    
                    
    # else :
    #     raise NotImplementedError
        
        
#Display of solution
# r_ways = route_ways[:10]


# ax=r_ways[0][0].plot()
# for route in r_ways:
#     route[0].plot(ax=ax)
#     route[-1].plot(ax=ax)
#     for point1, point2 in zip(route[:-1], route[1:]):
#         edges.LineSegment3D(point1, point2).plot(ax=ax, color = 'r')

# for list_jw in list_junc_wires :
#     for jw in list_jw :
#         ax = jw.path.plot(ax=ax)
        
    #     fix = point1
    #     # print(length, 'length')
    #     # print(points[0].point_distance(fix), points[-1].point_distance(fix))
        
    #     new_points = [points[0], fix, points[-1]]
    #     # abscissa_coord = [0, 0.5, 1]
    #     abscissa_coord = [0, fix_interval[0][0], 1]
        
    #     # for pt in new_points :
    #     #     pt.plot(ax=ax, color='m')
        
    #     tangents, points_jw = [], []
    #     for k in range(len(new_points)-1):
    #         p1, p2 = new_points[k:k+2]
    #         abs1, abs2 = abscissa_coord[k:k+2]
            
    #         tan1 = orls3d.point_at_abscissa((abs1+1/100)*length) - p1
    #         tan1.normalize()
            
    #         tan2 = orls3d.point_at_abscissa((abs2-1/100)*length) - p2
    #         tan2.normalize()
            
    #         tangents.append([tan1, tan2])
            
    #         points_jw.append([p1, p2])
        
    #     jw = []
    #     for pts, tans in zip(points_jw, tangents):
    #         lmin = pts[0].point_distance(pts[1])
    #         lmax= 3*lmin
        
    #         # ax=points[0].plot()
    #         # for pt in points[1:]:
    #         #     pt.plot(ax=ax)
    #         # for point1, point2 in zip(points[:-1], points[1:]):
    #         #     edges.LineSegment3D(point1, point2).plot(ax=ax, color = 'r')
    #         # edges.LineSegment3D(pts[0], pts[0]+0.1*tans[0]).plot(ax=ax,color='y')
    #         # edges.LineSegment3D(pts[1], pts[1]+0.1*tans[1]).plot(ax=ax,color='b')
    #         try:
    #             jw_to_add = wires.JunctionWire.curvature_radius(point1 = pts[0], 
    #                                                         tangeancy1 = tans[0], 
    #                                                         point2 = pts[1],
    #                                                         tangeancy2 = tans[1], 
    #                                                         targeted_curv = targeted_curv, 
    #                                                         diameter = diam, 
    #                                                         length_min = lmin, 
    #                                                         length_max = lmax)
    #             valid = True
    #         except ValueError:
    #             valid = False
                
    #         while not valid:
    #             lmax += 0.0001
    #             try:
    #                 jw_to_add = wires.JunctionWire.curvature_radius(point1 = pts[0], 
    #                                                         tangeancy1 = tans[0], 
    #                                                         point2 = pts[1],
    #                                                         tangeancy2 = tans[1], 
    #                                                         targeted_curv = targeted_curv, 
    #                                                         diameter = diam, 
    #                                                         length_min = lmin, 
    #                                                         length_max = lmax)
    #                 valid = True
    #             except ValueError:
    #                 valid = False
    #         # jw_to_add = wires.JunctionWire.curvature_radius(point1 = pts[0], 
    #         #                                                 tangeancy1 = tans[0], 
    #         #                                                 point2 = pts[1],
    #         #                                                 tangeancy2 = tans[1], 
    #         #                                                 targeted_curv = targeted_curv, 
    #         #                                                 diameter = diam, 
    #         #                                                 length_min = lmin, 
    #         #                                                 length_max = 3*lmin)
            
    #         # jw_to_add.path.plot(ax=ax)
            
    #         jw.append(jw_to_add)
              
        
    # else :
    #     print('length', orls3d.length())
    #     NotImplementedError
        
    # for junction in jw:
    #     junction.path.plot(ax=ax, color='k')
    #     print('radius_jw',junction.path.minimum_radius())


# # p1, p2 = route_ways[-2][0], route_ways[-2][-1]
# # t1 =  route_ways[-2][1] - route_ways[-2][0]
# # t1.normalize()
# # t2 = route_ways[-2][-2] - route_ways[-2][-1]
# # t1.normalize()
# # lmin = 50*p1.point_distance(p2)
# # lmax = 10*lmin
# # # try:
# # #     junction_wire = wires.JunctionWire.curvature_radius(p1, t1, p2, t2, 3*new_diameters[-2], new_diameters[-2], lmin, lmax)
# # #     valid = True
# # # except ValueError:
# # #     valid = False

# # # while not valid:
# # #     lmax += lmax*0.1
# # #     try:
# # #         junction_wire = wires.JunctionWire.curvature_radius(p1, t1, p2, t2, 3*new_diameters[-2], new_diameters[-2], lmin, lmax)
# # #         valid = True
# # #     except ValueError:
# # #         valid = False

# # junction_wire = wires.JunctionWire.curvature_radius(p1, t1, p2, t2, 3*new_diameters[-2], new_diameters[-2], lmin, lmax)
# # points =  route_ways[-2]
# # if len(points) == 3:
# #     bc = edges.BezierCurve3D(2, points)
# # else:
# #     bc = edges.BezierCurve3D(3, points)

# # ax=bc.plot(color = 'b')
# # for point1, point2 in zip(points[:-1], points[1:]):
# #     edges.LineSegment3D(point1, point2).plot(ax=ax, color = 'r')
# # for point in points:
# #     point.plot(ax=ax)
# # junction_wire.plot(ax=ax)
