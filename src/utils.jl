function save_sigma_profile(profile::SigmaProfile, filename::String)

    open(filename, "w") do io
        for i in eachindex(profile.sigma_grid)
            println(io,
                profile.sigma_grid[i],
                " ",
                profile.area_distribution[i])
        end
    end
end