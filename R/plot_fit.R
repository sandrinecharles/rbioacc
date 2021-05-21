

# df = Male_Gammarus_Single
# 
# load("data/fitMGS.rda")
# 
# fit
# 
# CGpred[t, 1] = (C0 - R[t]) * exp(-(E + M) * tp[t]) + R[t]
# 
# CGpred[t, 1] = (C0 - R[t] * (1 - exp((E + M)*tacc))) * exp(-(E + M) * tp[t])
# 
# ggplot(data = df) + 
#   theme_classic() +
#   geom_point(aes(x = time, y = conc))


