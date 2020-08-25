# cm inch transfer for matplotlib
def cm2inch(*tupl):
    inch = 2.54
    return tuple(i/inch for i in tupl)

# calculate the subplot size and plot height
# shape: two dimensional list, nrow and ncol
# margin: four dimensional list, left, bottom, right, top
# space: two dimensional list, height and width
def config_uniform_subplots( 
    plot_width = 9, 
    shape = [1, 1], 
    ratio = 0.8, 
    margin = [1.5, 1.2, 0.3, 0.3], 
    space = [1, 1] ): 

    subplot_width = ( plot_width 
                    - margin[0] - margin[2] 
                    - ( shape[1] - 1 ) * space[1] 
                    ) / shape[1]

    subplot_height = subplot_width * ratio

    plot_height = ( shape[0] * subplot_height
                  + margin[1] + margin[3]
                  + ( shape[0] -1 ) * space[0] )

    return plot_height, subplot_height, subplot_width
