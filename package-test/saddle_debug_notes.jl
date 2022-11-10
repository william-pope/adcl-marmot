# saddle_debug_notes.jl

# ISSUE: best value is to stay stopped at x=(10.984, 4.004, 1.738, 0.0)
#   - points to problem with solving/optimization
#   - approximate policy is able to find actions, so somehow no issue with neighboring nodes
#       - check their Q-values
#   - possible that forward neighbors are unsafe, while rear neighbors are fine?
#   - potentially related to dip in value function over velocity range

#   - nearest neighbor: ind_s_NN = 20523, x_NN = (11.0, 4.0, 1.5708, 0.0)
#   - see if any abnormalities in q_value_array at each neighbor
#       - 20522 -> 

#   - value at any v=0 node should be optimal reward-to-go, which happens when it speeds up
#   - q_value for v=0 actions should be q_value of best +0.5 action, minus one time step penalty

#   - seems like issue is due to heading angle?

# QUESTIONS
#   - is there a way to detect saddle points and keep it from getting stuck?
#       -

#   - can saddle points occur on all axes, or is this specific to the heading angle?
#       - actually a function of all states, because if vehicle was in a slightly diffrent spot, wouldn't be split on steering angles

#   - does this only occur when the vehicle is stopped?
#       - if vehicle was moving, 

#   - can HJB be solved with nearest neighbor approximations for the value?
#       - sort of works
#       - weird paths because for a lot of states, all steering angles land in the same grid
#           - this gives them all the same Q-value, so agent takes [-0.475, +0.5] whenever it can
#       - might work well for smaller grids, but that becomes intractable to solve
#       - takeaway: don't pursue NN solving any further

#   - at right heading angle, can get good V(s') by going right around obstacle
#   - at left heading angle, can get good V(s') by going left around obstacle
#   - in between heading angles, get somewhat worse V(s') because s' is close to obstacle
#       - however when taking Dv=0.0, s'=s
#       - interpolated V(s) is pretty good, because it takes the values from the left/right heading angles (which are truly good)
#   - if in-between angle was a grid node during solving, it would have a true worse value than left/right, which reflects that 
#       the only actions from that state are to go forward towards (or even in) the osbtacle

#   - in general, any state can get artifically inflated during interpolation
#       - there can be times where vehicle is moving (v!=0.0) and the state it lands in is worse than its current state
#       - issue arises when agent is able to stay in the same state
#       - think as long as vehicle keeps changing states, can't get stuck in artificial minima

#   - can vehicle  get itself into a minima? or does it only occur when the vehicle starts there?
#       - probably not in pure planning, because vehicle will never bring itself to a stop
#       - probably not in reactive planning, because RC will choose +0.5 eventually

#   - when vehicle is driving at v=0.5, will it ever choose Dv=-0.5 because minima state looks attractive?
#       - could be possible, but nearby states should be higher value -> chooses actions that leads to those
#       - even if lands there, should be able to push vehicle out
#       - agent shouoldn't be attracted to trough, but they are fairly wide regions so plausible that randomly sampling x_0 will produce states in troughs

#   - broadly, this is due to the value function being discontinuous (due to vehicle dynamics)

# FIXES
#   - gets stuck when staying stopped at current state has an artifically higher value than taking a forward action
#   - current state gets artifically inflated during interpolationn from neighboring states

#=
neighboring values =
    381.3764676899291
    966.7209427493063 (NN_1)
    -9202.264105107944
    -421.6617851131773
    975.2417002223822
    959.4486737530963 (NN_2)
    394.1617763642198
    -311.0306416955339
=#

# x=(10.984, 4.004, 1.738, 0.0)

#=
ind_s_NN_1 = 20523
x_NN_1 = (11.0, 4.0, 1.5708, 0.0)
Q-values = 
    -1.00000175e6
    -1.00000175e6
    -1.00000175e6
    -1.00000175e6
    -1.00000175e6
    -1.00000175e6
    -1.00000175e6
    964.7209427493063
    964.7209427493063
    964.7209427493063
    964.7209427493063
    964.7209427493063
    964.7209427493063
    964.7209427493063
    966.7209427493063 (*) -> [-0.475, +0.5], turning right and accelerating
    962.7298829869892
    959.0453341455855
    955.4736387001616
    941.40787648016
    930.6037309853538
    923.4259484595505
 =#

 # x=(10.984, 4.004, 1.738, 0.0)

#=
ind_s_NN_2 = 22204
x_NN_2 = (11.0, 4.0, 1.963, 0.0)
Q-values = 
    -1.00000175e6
    -1.00000175e6
    -1.00000175e6
    -1.00000175e6
    -1.00000175e6
    -1.00000175e6
    -1.00000175e6
    957.4486737530963
    957.4486737530963
    957.4486737530963
    957.4486737530963
    957.4486737530963
    957.4486737530963
    957.4486737530963
    878.6198292749152
    898.7484278557196
    920.4405117574911
    943.9803005725978
    948.9936787622885
    954.0666587504851
    959.4486737530963 (*) -> [+0.475, +0.5], turning left and accelerating
 =#

#=
generating Q-values between those heading angles

x_NN_1 = (11.0, 4.0, 1.5708, 0.0)
x_NN_2 = (11.0, 4.0, 1.963, 0.0)

x=11.0, y=4.0
theta | Q(s, a_15) | Q(s, a_21) | Q(s, stop)

1.5708 | 966.72 | 923.42 | 964.72
1.6100 | 964.03 | 916.46 | 963.99   (near saddle)

1.6500 | 960.37 | 913.99 | 963.25   (saddle)
1.7500 | 901.54 | 927.60 | 961.40   (saddle)
1.8500 | 874.10 | 947.86 | 959.55   (saddle)

1.9500 | 876.20 | 958.17 | 957.70   (near saddle)
1.9630 | 878.52 | 959.40 | 957.46