#' Create a custom legend
#'
#' A customized legend sections the group taxonomy as headers
#' and lists the subgroup taxonomies below each group header, for all the groups.
#'
#' @param mdf data.frame melted dataframe
#' @param cdf data.frame color dataframe
#' @param group_level string of larger taxonomic group
#' @param subgroup_level string of smaller taxonomic group
#' @param x x-value plotted
#' @param y y-value plotted
#' @param legend_key_size numeric determines overall size of the legend keys
#' @param legend_text_size integer determines size of legend text
#' @param legend_orientation string determines vertical or horizontal orientation of legend
#'
#' @import dplyr
#' @import ggplot2
#' @import cowplot
#'
#' @return complete_legend ggplot
#' @export
#'
#' @examples
#' legend_new <- custom_legend(mdf_GP, cdf_GP)
#'
#'
custom_legend <- function (mdf, cdf, group_level = "Phylum", subgroup_level = "Genus", x = "Sample",
                           y = "Abundance", legend_key_size = 0.4, legend_text_size = 10, legend_orientation = "vertical")
{
    if (is.null(mdf[[group_level]])) {
        stop("mdf 'group_level' does not exist")
    }
    
    if (is.null(mdf[[subgroup_level]])) {
        stop("mdf 'subgroup_level' does not exist")
    }
    
    if (is.null(cdf$hex)) {
        stop("cdf 'hex' does not exist")
    }
    
    if (!(legend_orientation %in% c("vertical", "horizontal"))) {
        stop("legend orientation must be \"vertical\" or \"horizontal\"")
    }
    
    col_name_group <- paste0("Top_", group_level)
    col_name_subgroup <- paste0("Top_", subgroup_level)
    
    group_level_names <- unique(cdf[[col_name_group]])
    
    for (i in 1:length(group_level_names))
    {
        if( i == 1)
        {
            complete_legend <-individual_legend (mdf = mdf, cdf = cdf, group_name = group_level_names[i], col_name_group = col_name_group, col_name_subgroup = col_name_subgroup, x = x, y = y, legend_key_size = legend_key_size, legend_text_size = legend_text_size)
        }
        else
        {
            new_legend <-individual_legend (mdf = mdf, cdf = cdf, group_name = group_level_names[i], col_name_group = col_name_group, col_name_subgroup = col_name_subgroup, x = x, y = y, legend_key_size = legend_key_size, legend_text_size = legend_text_size)
            
            complete_size <- i -1
            new_size <- 1
            
            if (legend_orientation == "vertical") {
                complete_legend <-plot_grid(complete_legend, new_legend, ncol = 1, rel_heights = c(complete_size,new_size))
            } else if (legend_orientation == "horizontal") {
                complete_legend <-plot_grid(complete_legend, new_legend, ncol = 2, rel_widths = c(complete_size,new_size))
            }
        }
    }
    complete_legend
}


#' Creates a small custom legend for one group variable
#'
#' Create a custom legend for one individual group, with the subgroups below
#'
#' @param mdf data.frame melted dataframe
#' @param cdf data.frame color dataframe
#' @param group_name name of the larger taxonomic group to extract small legend from
#' @param col_name_group string column name of larger taxonomic group
#' @param col_name_subgroup string column nameof smaller taxonomic group
#' @param x x-value plotted
#' @param y y-value plotted
#' @param legend_key_size numeric determines overall size of the legend keys
#' @param legend_text_size integer determines size of legend text
#'
#' @import dplyr
#' @import ggplot2
#' @import cowplot
#'
#' @return legend ggplot
#' @keywords internal
#'
#'
#'
individual_legend <- function (mdf,
                               cdf,
                               group_name,
                               col_name_group = "Top_Phylum",
                               col_name_subgroup = "Top_Genus",
                               x = "Sample",
                               y = "Abundance",
                               legend_key_size = 0.4,
                               legend_text_size = 10)
{
    select_mdf <- mdf %>% filter(!!sym(col_name_group) == group_name)
    select_cdf <- cdf %>% filter(!!sym(col_name_group) == group_name)
    
    select_plot <- ggplot(select_mdf,
                          aes(x = .data[[x]], y = .data[[y]], color = .data[[col_name_subgroup]], text = .data[[col_name_subgroup]])) +
        geom_point() +
        scale_color_manual(name = group_name,
                          values = select_cdf$hex,
                          breaks = select_cdf[[col_name_subgroup]]) +
        theme(legend.justification = "left") +
        theme(legend.title = element_text(face = "bold")) +
        theme(legend.key.size = unit(legend_key_size, "lines"), text=element_text(size=legend_text_size))
    
    legend <- cowplot::get_legend(select_plot)
}