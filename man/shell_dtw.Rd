\name{shell_dtw}
\alias{shell_dtw}
\title{Dynamic time Warping of Shell records}
\usage{
shell_dtw(D18o,MgCa)
}
\description{
Produces correlation graph of both proxies and some basic statistics.
It is important that the two records are split into two dataframes or tibbles with specific column names.
The D18o dataframe has one column named `sample_seq`, which contains the sample sequence as integer numbers,
and one column called `d18o`, which contains the corresponding d18o values.
The MgCa dataframe has one column named `dist`, which contains the distance of the samples to the shell edge and one column called `mg_ca`, which contains the corresponding Mg/Ca values.
(if the distances are not included, they can be calculated via:
  `MgCa %>%
  mutate(dist=row_number()*0.03))`

}
\examples{

MgCa <- read_csv("FR-P_20210515_G30_ANT_LINE-SCAN_analyzed_data.csv", col_names = FALSE) %>%
  rename(x=X1,y=X2,z=X3,mg_ca=X4,std=X5, rel_std=X6) %>%
  mutate(dist=row_number()*0.03) # adds the distance values if samples are 30 µm (0.03 mm) apart

D18o <- tibble::tribble(
          ~sample_seq,  ~d18o,
           1L, 0.85,
           2L, 1.03,
           3L, 1.03,
           4L, 0.91,
           5L, 0.78,
           6L, 0.83,
           7L, 1.23,
           8L, 1.29,
           9L, 1.39,
          10L, 1.54,
          11L, 1.58,
          12L, 1.39,
          13L, 1.43,
          14L, 1.08,
          15L, 0.99,
          16L, 1.01,
          17L, 0.81,
          18L, 0.76,
          19L, 0.52
          )

shell_dtw(D18o,MgCa)

}
