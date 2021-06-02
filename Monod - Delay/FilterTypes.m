classdef FilterTypes < Simulink.IntEnumType
  enumeration
    EKF(0)
    iEKF(1)
    UKF(2)
    PF(3)
  end
end