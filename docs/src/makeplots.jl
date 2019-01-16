using Plots, LaTeXStrings
Plots.gr()

function  example_plot_function(plot_dir)
    # Simulate Keplerian Orbit
    plot_name = "example_plot"
    println("Generating: $plot_name.svg")

    # Main plotting routine
    x = 1:10; y = rand(10,2) # 2 columns means two lines
    plot(x,y,title="Two Lines",label=["Line 1" "Line 2"],lw=3)

    # Save figure
    Plots.savefig(plot_dir * "/$plot_name.svg")
end

# Generate all plots
function makeplots()
    # Start building plots
    println("Generating plots")

    # Construct output directory path
    plot_dir = (pwd()[end-3:end] == "docs") ? "build/plots" : "docs/build/plots"

    # Remove directory if it already exists to ensure clean build
    if isdir(plot_dir)
        rm(plot_dir, recursive=true)
    end

    # Create output directory
    mkdir(plot_dir)

    # Call plotting functions
    example_plot_function(plot_dir)
end