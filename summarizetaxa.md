Creating a summarize taxa function
================

# Open necessary libraries

This code uses the phyloseq S4 class object created in the phyloseq
pipeline (see: <https://benjjneb.github.io/dada2/tutorial.html>) to
summarize taxa in a function similar to how QIIME does. It also uses
data.table to quickly make data into a table (this was more a past
concern when the code was originally written, now dplyr works just as
fast as data.table in most instances). This code was first mentioned in
comment: <https://github.com/joey711/phyloseq/issues/818> but the
GroupBy function was not working as intended. This solves as my bandaid
fix to trying to figure out what was going on with this issue.

``` r
library("phyloseq")
library("data.table")
```

# Create the fast\_melt function

The `fast_melt` function will be used to coerce the OTU table with raw
counts and the tax table with taxa names based on OTUs from the phyloseq
object to be able to create a data table with counts for
taxa.

``` r
fast_melt = function(physeq){ # converts raw otu counts to taxa with values
  # supports "naked" otu_table as `physeq` input.
  otutab = as(otu_table(physeq), "matrix")
  if(!taxa_are_rows(physeq)){otutab <- t(otutab)} # imports OTU names based on the rows
  otudt = data.table(otutab, keep.rownames = TRUE) # creates a data table of the otus, thought that data.table will be faster than dplyr
  setnames(otudt, "rn", "taxaID") # allows for taxa IDs to be associated with the otu number
  # Enforce character taxaID key
  otudt[, taxaIDchar := as.character(taxaID)] # looks at the column to name the taxa
  otudt[, taxaID := NULL]
  setnames(otudt, "taxaIDchar", "taxaID")
  # Melt count table
  mdt = melt.data.table(otudt, 
                        id.vars = "taxaID",
                        variable.name = "SampleID",
                        value.name = "count")
  # Remove zeroes, does NA calculations for reporting
  mdt <- mdt[count > 0][!is.na(count)] 
  # Calculate relative abundance
  mdt[, RelativeAbundance := count / sum(count), by = SampleID] # creates a data table out of 100
  if(!is.null(tax_table(physeq, errorIfNULL = FALSE))){
    # If there is a tax_table, join with it. Otherwise, skip this join.
    taxdt = data.table(as(tax_table(physeq, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE) # identifies the taxa table and imports it as a matrix from your phyloseq data
    setnames(taxdt, "rn", "taxaID") # setnames of the objects that we want to call on later in the function
    # Enforce character taxaID key
    taxdt[, taxaIDchar := as.character(taxaID)] # the := operator is used in data.table to identify a new column
    taxdt[, taxaID := NULL] # creates a blank column for the taxaID 
    setnames(taxdt, "taxaIDchar", "taxaID") # This is a convenience function that sets the names on an object and returns the object. It is most useful at the end of a function definition where one is creating the object to be returned and would prefer not to store it under a name just so the names can be assigned (in this case, taxaIDs and taxaIDchar)
    # Join with tax table
    setkey(taxdt, "taxaID") # setkey sorts a data.table and marks it as sorted with an attribute
    setkey(mdt, "taxaID") # setkey used to keep the taxaID with relative abundance
    mdt <- taxdt[mdt] # merging the taxaIDs with relative abundance
  }
  return(mdt) # returns the data table that contains relative abundance and taxaIDs
}
```

# Create summarize\_taxa function

The purpose of this is taking the phyloseq object and creating a summary
of taxa. Previously, this did not correctly have a GroupBy function. It
would still add up to 1 even if the samples were supposed to be broken
into two groups (i.e., if you had 1 factor containing 2 levels that you
wanted to group, it would divide everything by 2 to still add everything
up to 1). This was counterproductive to what we wanted the GroupBy
function to work, so we created a statement for if the GroupBy was added
to correctly add up.

``` r
# Actual summarize taxa function
summarize_taxa = function(physeq, Rank, GroupBy = NULL){ # where you can put in the phyloseq object, the rank of interest, and then groupby for a variable
  Rank <- Rank[1] # pick a rank of interest 
  if(!Rank %in% rank_names(physeq)){ # based on rank names available from your phyloseq 
    message("The argument to `Rank` was:\n", Rank, # error message if you type your name wrong
            "\nBut it was not found among taxonomic ranks:\n",
            paste0(rank_names(physeq), collapse = ", "), "\n", # allows you to see the ranks you can choose from
            "Please check the list shown above and try again.")
  }
  if(!is.null(GroupBy)){ # pick a variable to group by
    GroupBy <- GroupBy[1] # Must be a character
    if(!GroupBy %in% sample_variables(physeq)){ # Use sample variables from phyloseq object
      message("The argument to `GroupBy` was:\n", GroupBy, # gives you the error message if your groupby doesn't exist and provides a list of options
              "\nBut it was not found among sample variables:\n",
              paste0(sample_variables(physeq), collapse = ", "), "\n",
              "Please check the list shown above and try again.")
    }
  }
  # Start with fast melt function first defined
  mdt = fast_melt(physeq) # This takes your phyloseq OTU and tax tables and 
  if(!is.null(GroupBy)){
    # Add the variable indicated in `GroupBy`, if provided.
    sdt = data.table(SampleID = sample_names(physeq), # SampleID is your library indicator
                     var1 = get_variable(physeq, GroupBy)) # imports your variable 
    setnames(sdt, "var1", GroupBy) # creates columns that we will put data into later
    # Join using the setkey function as defined in the data.table package
    setkey(sdt, SampleID) # sampleID for the variables in the phyloseq object
    setkey(mdt, SampleID) # sampleID for the melted phyloseq object containing transformed otu counts with taxa names
    mdt <- sdt[mdt] # merge into the data table
  }
  # Summarize
  if(!is.null(GroupBy)){ 
    var1 = get_variable(physeq, GroupBy) # retrieves the variable from the phyloseq object
    nvars = nlevels(var1) # counts how many variables there are 
    Nsamples = nsamples(physeq) # Phyloseq sample counts
    summarydt = mdt[, list(meanRA = (sum(RelativeAbundance)/ Nsamples) * nvars, # Multiplies by the number of levels for the specific variable
                           sdRA = sd(RelativeAbundance),
                           minRA = min(RelativeAbundance),
                           maxRA = max(RelativeAbundance)),
                    by = c(Rank, GroupBy)]
  } else {
    Nsamples = nsamples(physeq) # Nsamples is a phyloseq function that will pull the number of samples from phyloseq
    # This is used to remove the issue of samples with 0's that will falsely skew mean relative abundance
    summarydt = mdt[, list(meanRA = sum(RelativeAbundance)/Nsamples, # calculate the meanRA based on the number of samples
                           sdRA = sd(RelativeAbundance), # identifies the standard deviation
                           minRA = min(RelativeAbundance), # minimum 
                           maxRA = max(RelativeAbundance)), # maximum
                    by = c(Rank, GroupBy)]
  return(summarydt) # returns the summarized data table
}
}
```

Some things to consider are that the min and max RA have been previously
tested and shown to be overly variable in comparison to the true data
(on the to-do list to fix), however the meanRAs have been tested against
QIIME2 outputs to confirm that you would get similar values to
summarizetaxa.py. The purpose of this function is to be able to report
relative abundance for manuscript writing.
