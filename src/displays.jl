WIDTH = 67

function display_head()
    println(repeat("*",WIDTH))
    println("B&B-screening solver for the L0-penalized least-square problem")
    println("Version : $(version())")
    println("Author  : $(author())")
    println("Contact : $(contact())")
    println("License : $(license())")
    println()
    println(repeat("-",WIDTH))
    println("|        Node         |        Bounds       |        Branch       |")
    println("|  nexpl  depth  queue|     lb     ub    gap|     S0     S1   Sbar|")
    println(repeat("-",WIDTH))
    return nothing
end
    
function display_node(tree::Tree, node::Node)

    @printf(
        "|%7d%7d%7d|%7.3f%7.3f%7.0e|%7d%7d%7d|",
        tree.nexpl,
        depth(node),
        length(tree.queue),
        best_lb(tree),
        tree.ub,
        global_gap(tree),
        length(node.branch.S0),
        length(node.branch.S1),
        length(node.branch.Sbar),
    )
    println()
    return nothing
end

function diplay_results(results::BnbResults)
    println("B&B results")
    println("  Status     : $(results.termination_status)")
    println("  Solve time : $(results.solve_time)")
    println("  Nodes      : $(results.node_count)")
    println("  Objective  : $(results.objective_value)")
    return nothing
end

function display_tail(results::BnbResults)
    println(repeat("-",WIDTH))
    println()
    diplay_results(results)
    println()
    println(repeat("*",WIDTH))
    return nothing
end
