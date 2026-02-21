# placeholder for if in the future I'm to implement the full COSMOspace EQUATIONs 

struct MixtureComponent
    profile::SigmaProfile
    mole_fraction::Float64
end

struct Mixture
    components::Vector{MixtureComponent}
end