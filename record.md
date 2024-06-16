先前：无滑移边界条件，得不到好的结果
Roe1：反射边界条件，收敛，Pressure Ratio before and after shockwave:  2.106214421687035
Roe2：前段固壁u=2,v=0，后段反射边界条件，收敛，Pressure Ratio before and after shockwave:  2.1060402617365552
Roe3：发现了一个错误，修正后，收敛，Pressure Ratio before and after shockwave： 2.1249888558074908
Roe4: 调参Pressure Ratio before and after shockwave:  2.1265850998796445
Roe5：继续调参（更改熵修正方式）：Pressure Ratio before and after shockwave:  2.1261944235859183，耗散大

采用Roe4稍作修改。
Pressure Ratio before and after shockwave:  2.1276712988414763
Mach Number after shockwave:  1.4579336952520368

Roe6：修改网格，激波很sharp
Pressure Ratio before and after shockwave:  2.1630171151043713
Mach Number after shockwave:  1.4495422724299725

Roe7：修改时间步长
Pressure Ratio before and after shockwave:  2.1634308224109624
Mach Number after shockwave:  1.4500208352507493



