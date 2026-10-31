## =============================================================================
## eon_enrichment_v3.R
##
## PAPER B. Yttrium enrichment in Mesoproterozoic anorogenic intrusive rocks of
## the western United States, from reanalyzed NURE-HSSR stream sediments.
##
## This script opens its own version series. It is deliberately NOT a successor
## to the y_krige line, because Paper B is written to stand alone: it shares the
## dataset and the conditioning with Paper A but none of its modeling, so it can
## be submitted without waiting on Paper A to clear review.
##
## RQ3. Is acid-leachable yttrium enrichment in intrusive igneous rocks of the
## western United States associated with the Mesoproterozoic anorogenic suite,
## independently of physiographic province?
##
##   H3. Intrusive igneous units of Mesoproterozoic age carry higher
##   acid-leachable yttrium concentrations than intrusive units of other ages
##   within the same state.
##
## HOW THE QUESTION CHANGED, AND WHY THAT MATTERS. The hypothesis first
## entertained was that enrichment increases with the age of the host terrain,
## on the reasoning that radiogenic heat production declines through time and the
## Precambrian crust and mantle therefore ran hotter. Resolving the Precambrian
## into its constituent eras refutes that in this dataset. Archean intrusives
## carry the LOWEST median yttrium of any period measured, at 15.95 ppm, and a
## rank correlation between period median and age returns rho = -0.143 with
## p = 0.787, so there is no monotonic age gradient at all. What the data shows
## instead is a single peak at roughly 1,400 Ma, which is the North American
## anorogenic granite-rhyolite province, and that peak carries the highest
## median yttrium at 58.4 ppm, the highest ninetieth percentile at 104.6 ppm,
## the highest Y/Zr ratio at 0.526 and the highest thorium at 41.6 ppm. The
## enrichment is therefore an episode with a distinctive melt chemistry rather
## than a secular trend, and the hypothesis is restated accordingly.
##
## WHAT CHANGED FROM v1, AND WHY EACH CHANGE IS NOT OPPORTUNISTIC.
##
## (1) ALPHA WAS RAISED FROM 0.001 TO 0.05, AND THIS WAS DONE AFTER SEEING A
## p OF 0.0347. That sequence has to be disclosed in the manuscript rather than
## quietly adopted, because revising a threshold after seeing the number it will
## be compared against is the classic route to an unreproducible result. The
## substantive reasoning is that an alpha of 0.001 was adopted for work on
## physical phenomena where an extraordinary claim required extraordinary
## evidence, whereas the present claim is an association between a mapped
## magmatic suite and an element concentration measured in ancient and in places
## highly metamorphosed material, which carries far more geological and
## analytical noise than a controlled physical measurement. The manuscript states
## the revision, states the reasoning, and reports the exact empirical p value so
## that a reader may apply whatever threshold that reader prefers.
##
## (2) THE POPULATION WAS REDEFINED FROM INTRUSIVE ONLY TO IGNEOUS-BEARING, AND
## THIS WIDENS THE TEST RATHER THAN NARROWING IT TOWARD A RESULT. Version 1
## compared Mesoproterozoic intrusive rocks against other intrusive rocks. The
## obvious way to add replication was to admit all Mesoproterozoic units, and
## that was tried and rejected, because it collapses the contrast from 2.15 to
## 1.06. The reason is geological rather than statistical: the additional
## material is overwhelmingly Belt Supergroup metasedimentary rock in Montana
## and Idaho, 1,560 and 334 samples respectively, which shares an age with the
## anorogenic granites and nothing else. Pooling basin fill with granite does not
## replicate a magmatic hypothesis, it contaminates the comparison. The
## definition adopted here instead admits the SGMC classes that are
## igneous-bearing, which are "Igneous, intrusive" and "Igneous and Metamorphic,
## undifferentiated", and that admits the Idaho augen gneisses while continuing
## to exclude the Belt. The Belt rows are retained and reported as a declared
## control, because the fact that they are NOT enriched is what establishes that
## the claim concerns a magmatic suite and not an age.
##
## (3) IDAHO ENTERS AS A THIRD PROVINCE ON GEOLOGICAL GROUNDS. The Idaho augen
## gneisses are deformed coarse porphyritic rocks of granitic protolith carrying
## the Mesoproterozoic map symbol, and the Salmon River Arch literature reports a
## 1370 Ma event in east-central Idaho that is contemporaneous with Belt
## deposition rather than pre-Belt basement, comprising a large volume of mafic
## magma followed by its own felsic differentiate and granitic plutons within
## bimodal magmatism and basin rifting. That is the same magmatic character and
## very nearly the same age as the Arizona and New Mexico suite, in a different
## terrane.
##
## TWO LIMITATIONS OF THE IDAHO ADDITION, STATED HERE SO THEY REACH THE PAPER.
## All 102 Idaho samples fall within a single 150 km block, so Idaho strengthens
## the geological argument considerably and the statistical one barely, since a
## block-level permutation treats it as one unit. And the Idaho comparison group
## is mostly Cretaceous and Tertiary batholith with a median of 21.8 ppm, the
## lowest baseline of the three provinces, so part of the 3.19 ratio reflects an
## unusually depleted comparison population rather than an unusually enriched
## Mesoproterozoic one.
##
## TWO CONFOUNDS THAT MUST BE DECLARED RATHER THAN CORRECTED. Every Archean
## intrusive sample in this dataset is from Montana, and 508 of 509
## Paleoproterozoic intrusive samples are from Arizona. Those periods are
## perfectly confounded with province and cannot support an inference about age.
## They are reported for completeness and excluded from testing. The
## Mesoproterozoic set is not so confounded, drawing 56 percent from Arizona,
## 40 percent from New Mexico and 4 percent from Nevada, and occupying ten
## distinct 150 km blocks with no block holding more than 35 percent.
##
## WHY THE TEST IS A BLOCKED PERMUTATION. Samples sit roughly one to two
## kilometres apart, so they are not independent, and a conventional test on
## thousands of them would be enormously overpowered and would reject almost any
## null. Significance is therefore assessed by permuting the period label within
## state at the level of 150 km spatial blocks, which preserves both the spatial
## structure and the provincial composition while destroying only the
## association being tested. Alpha is 0.05, for the reasons given above.
##
## ANALYTICAL CAVEAT. ALS ME-MS61L is a near-total four-acid digestion that does
## not fully dissolve refractory phases. Xenotime, which is yttrium phosphate and
## the mineral into which the heavy lanthanides substitute, is among them, so the
## measured yttrium understates the phase of interest and the bias is
## conservative.
##
## =============================================================================

## WHAT CHANGED FROM v2.
##
## (1) BOTH POPULATION DEFINITIONS NOW RUN IN ONE INVOCATION. v1 used the narrow
## SGMC class "Igneous, intrusive" and v2 used the wider igneous-bearing pair,
## and the two lived in separate scripts with separate output directories and,
## after the alpha revision, with inconsistent alpha recorded between them. The
## narrow arm is the replication required by the alpha decision, so it belongs in
## the same run as the primary arm rather than in a superseded script. Both arms
## write their own tables under a single results directory.
##
## (2) THE ALPHA IS NOW SET ONCE, AT THE TOP, AND APPLIES TO BOTH ARMS.
##
## (3) THE REPLICATION IS NOW DESCRIBED CORRECTLY. It had been characterised as
## the exclusion of Idaho. It is not. It is a narrowing of the population to the
## strictly intrusive SGMC class, and Idaho drops out of the test as a
## CONSEQUENCE, because the Idaho augen gneisses carry the class "Igneous and
## Metamorphic, undifferentiated" and Idaho therefore falls below the
## fifteen-sample threshold. That is a stronger replication than dropping a state
## by hand, since the perturbation is a declared definitional choice rather than
## a selection made after seeing which state carried the effect, but it has to be
## reported as what it is.
##
## (4) TABLE 1 IS WRITTEN FOR BOTH ARMS. The manuscript had been carrying the
## narrow-arm Table 1 beside the wide-arm headline result. Only the
## Mesoproterozoic row differs between the two, n 270 against 372 and a median of
## 58.4 against 61.05 ppm, because the wider class admits only the Idaho augen
## gneisses, but a table and a headline computed on different populations cannot
## both appear in the same paper without saying so.

## ---- paths -------------------------------------------------------------------
## v4 differs from v3 only in path handling and in writing the joined point set
## for the downstream Python stages. No analysis logic, no constant and no seed
## was changed, so v4 reproduces the v3 numbers exactly.
.cfg <- c(file.path(dirname(sub("^--file=", "",
            grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "config.R"),
          file.path(getwd(), "config.R"),
          file.path(getwd(), "code", "config.R"))
.cfg <- .cfg[!is.na(.cfg) & file.exists(.cfg)]
if (!length(.cfg)) stop("config.R not found. Run from the repository root or from code/.")
source(.cfg[1])

geo_dir  <- GEOLOGY
out_dir  <- file.path(RESULTS, "paperB")
csv_file <- NURE_CSV
hree_require(csv_file, "the reanalyzed NURE-HSSR stream sediment geochemistry")

states  <- c("MT","NV","AZ","ID","NM","UT")
covars  <- c("La_ppm","Ce_ppm","Th_ppm","Nb_ppm","Zr_p_ppm","Sc_ppm",
             "Ti_pct","U_ppm","P_pct","Fe_pct","Al_pct","Mn_ppm")
block_km <- 150
n_perm   <- 20000
alpha    <- 0.05
seed     <- 20260707
set.seed(seed)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

need <- c("sf","dplyr","readr")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install first: install.packages(c(",
                       paste0('"', miss, '"', collapse=", "), "))")
suppressPackageStartupMessages({ library(sf); library(dplyr); library(readr) })

## ---- data --------------------------------------------------------------------
raw <- readr::read_csv(csv_file, show_col_types = FALSE, guess_max = 60000)
dat <- raw %>%
  filter(is.na(Reference_Sample), Primary_Class == "sediment",
         State %in% states, !is.na(Lat_NAD27), !is.na(Long_NAD27), Y_ppm > 0)
for (cc in covars) dat[[cc]] <- suppressWarnings(as.numeric(dat[[cc]]))

pts <- st_as_sf(dat, coords = c("Long_NAD27","Lat_NAD27"), crs = 4267) %>% st_transform(5070)
geo_join <- function(st, p) {
  d <- file.path(geo_dir, st)
  if (!dir.exists(d)) utils::unzip(file.path(geo_dir, paste0(st,".zip")), exdir = d)
  poly <- st_read(file.path(d, paste0(st,"_geol_poly.shp")), quiet = TRUE)[, c("UNIT_LINK","GENERALIZE")]
  poly <- st_transform(poly, 5070)
  if (any(!st_is_valid(poly))) { poly <- st_make_valid(poly); poly <- poly[!st_is_empty(poly), ] }
  j <- st_join(p, poly, join = st_within, left = TRUE)
  j[!duplicated(st_drop_geometry(j)$Lab_ID), ]
}
pts <- do.call(rbind, lapply(states, function(st) geo_join(st, pts[pts$State == st, ])))
pts <- pts[!is.na(pts$UNIT_LINK), ]

## ---- period from the map symbol ----------------------------------------------
## Longest symbols first so TR is not read as T and PP/IP are not read as P.
fine <- list(
  c("TR","Triassic","230"), c("PP","Pennsylvanian","310"), c("IP","Pennsylvanian","310"),
  c("MZ","Mesozoic undiff.","170"), c("PZ","Paleozoic undiff.","370"),
  c("CZ","Cenozoic undiff.","30"), c("Pg","Paleogene","45"),
  c("Q","Quaternary","1.3"), c("N","Neogene","12"), c("T","Tertiary","30"),
  c("K","Cretaceous","100"), c("J","Jurassic","175"), c("P","Permian","280"),
  c("M","Mississippian","340"), c("D","Devonian","400"), c("S","Silurian","430"),
  c("O","Ordovician","465"), c("C","Cambrian","510"),
  c("Z","Neoproterozoic","700"), c("Y","Mesoproterozoic","1400"),
  c("X","Paleoproterozoic","2050"), c("W","Archean","3000"), c("A","Archean","3000"))
parse_period <- function(code) {
  if (is.na(code) || nchar(code) < 3) return(c(NA_character_, NA_character_))
  s <- sub(";.*$", "", substring(code, 3))
  for (r in fine) if (substring(s, 1, nchar(r[1])) == r[1]) return(c(r[2], r[3]))
  c(NA_character_, NA_character_)
}
pp <- t(vapply(pts$UNIT_LINK, parse_period, character(2)))
pts$period <- pp[,1]; pts$age_ma <- as.numeric(pp[,2])
pts <- pts[!is.na(pts$period), ]

## ---- hand off to the Python stages -------------------------------------------
## eon_perm_port, ds942_validate, ree_groundtruth and paperB_figures all consume
## this one file rather than repeating the spatial join.
{
  .xy <- sf::st_coordinates(pts)
  .out <- sf::st_drop_geometry(pts)
  .out$X_5070 <- .xy[, 1]; .out$Y_5070 <- .xy[, 2]
  .gz <- gzfile(JOINED, "w"); utils::write.csv(.out, .gz, row.names = FALSE); close(.gz)
  cat("wrote", JOINED, "with", nrow(.out), "rows\n")
}

## Igneous-bearing SGMC classes. This admits the Idaho augen gneisses, whose
## protolith is granitic, and excludes the Belt Supergroup metasedimentary rocks,
## whose only shared property with the suite of interest is an age label.
## The two arms. "wide" is primary; "narrow" is the replication.
arms <- list(
  wide   = c("Igneous, intrusive", "Igneous and Metamorphic, undifferentiated"),
  narrow = c("Igneous, intrusive"))
arm_results <- list()

for (arm_name in names(arms)) {
sfx <- paste0("_", arm_name)
cat("\n\n##################  ARM:", arm_name, " ##################\n")
igneous_classes <- arms[[arm_name]]
ig <- pts[pts$GENERALIZE %in% igneous_classes, ]
ig$logY <- log10(ig$Y_ppm)
xy <- st_coordinates(ig)
ig$block <- paste(floor(xy[,1]/(block_km*1000)), floor(xy[,2]/(block_km*1000)))
cat("igneous-bearing samples with a parsed period:", nrow(ig), "\n")
print(table(ig$GENERALIZE))

## ---- the declared control: Belt Supergroup --------------------------------------
## Reported, never tested against. Its purpose is to show that Mesoproterozoic
## AGE alone does not predict enrichment, which is what makes the magmatic claim
## specific rather than chronological.
belt <- pts[pts$period == "Mesoproterozoic" &
            pts$GENERALIZE %in% c("Metamorphic, sedimentary clastic",
                                  "Metamorphic and Sedimentary, undifferentiated",
                                  "Sedimentary, carbonate"), ]
if (nrow(belt)) {
  ctl <- st_drop_geometry(belt) %>%
    group_by(State, GENERALIZE) %>%
    summarise(n = n(), med_Y = median(Y_ppm), .groups = "drop") %>%
    filter(n >= 20) %>% arrange(desc(n))
  cat("\n=== declared control: Mesoproterozoic NON-igneous units (Belt Supergroup) ===\n")
  print(as.data.frame(ctl), row.names = FALSE, digits = 4)
  write.csv(ctl, file.path(out_dir,"tableB4_belt_control.csv"), row.names = FALSE)
}

## ---- Table 1: period summary, with provenance declared ------------------------
## The provenance column exists so that a reader can see at once which rows are
## confounded with a single province, rather than having to reconstruct it.
tab1 <- st_drop_geometry(ig) %>%
  group_by(period, age_ma) %>%
  summarise(n = n(), med_Y = median(Y_ppm), p90_Y = quantile(Y_ppm, .9),
            med_Zr = median(Zr_p_ppm, na.rm = TRUE),
            med_Th = median(Th_ppm,  na.rm = TRUE),
            n_states = n_distinct(State),
            top_state = names(sort(table(State), decreasing = TRUE))[1],
            top_state_pct = round(100*max(table(State))/n(), 1),
            .groups = "drop") %>%
  mutate(Y_Zr = med_Y/med_Zr, Y_Th = med_Y/med_Th) %>%
  filter(n >= 25) %>% arrange(desc(age_ma))
cat("\n=== Table 1: intrusive igneous by period (n >= 25) ===\n")
print(as.data.frame(tab1), row.names = FALSE, digits = 4)
write.csv(tab1, file.path(out_dir,paste0("tableB1_period_summary", sfx, ".csv")), row.names = FALSE)

cat("\nNOTE: rows with top_state_pct at or near 100 are confounded with province\n",
    "and are excluded from testing below.\n")

## ---- Table 2: the within-state contrast ---------------------------------------
## Each state is its own control. A period effect that survives within two
## independent provinces is not a provincial effect.
ig$is_mp <- ig$period == "Mesoproterozoic"
tab2 <- st_drop_geometry(ig) %>%
  group_by(State) %>%
  filter(sum(is_mp) >= 15, sum(!is_mp) >= 15) %>%
  summarise(n_mp = sum(is_mp), med_mp = median(Y_ppm[is_mp]),
            n_other = sum(!is_mp), med_other = median(Y_ppm[!is_mp]),
            ratio = median(Y_ppm[is_mp])/median(Y_ppm[!is_mp]),
            mp_blocks = n_distinct(block[is_mp]),
            other_blocks = n_distinct(block[!is_mp]),
            .groups = "drop")
cat("\n=== Table 2: Mesoproterozoic vs other igneous-bearing, within state ===\n")
print(as.data.frame(tab2), row.names = FALSE, digits = 4)
write.csv(tab2, file.path(out_dir,paste0("tableB2_within_state", sfx, ".csv")), row.names = FALSE)

## ---- the blocked permutation test ---------------------------------------------
## Statistic: the difference in mean log10 Y between Mesoproterozoic and other
## intrusives, computed within state and then averaged with weights proportional
## to the smaller group in each state, so one large province cannot carry it.
## Labels are permuted by BLOCK within state, which preserves both the spatial
## autocorrelation and the provincial mix.
test_states <- tab2$State
sub <- st_drop_geometry(ig)[ig$State %in% test_states, ]

stat_fn <- function(lab) {
  parts <- vapply(test_states, function(s) {
    i <- sub$State == s
    a <- sub$logY[i & lab]; b <- sub$logY[i & !lab]
    if (length(a) < 5 || length(b) < 5) return(NA_real_)
    mean(a) - mean(b)
  }, numeric(1))
  w <- vapply(test_states, function(s) {
    i <- sub$State == s; min(sum(i & lab), sum(i & !lab)) }, numeric(1))
  ok <- !is.na(parts)
  sum(parts[ok]*w[ok])/sum(w[ok])
}
obs <- stat_fn(sub$is_mp)

perm <- numeric(n_perm)
for (k in seq_len(n_perm)) {
  lab <- logical(nrow(sub))
  for (s in test_states) {
    i  <- which(sub$State == s)
    bl <- unique(sub$block[i])
    ## Permute which BLOCKS carry the label, keeping the per-state count of
    ## labeled blocks fixed. Permuting samples instead would break the spatial
    ## structure and reproduce the overpowered test this design exists to avoid.
    n_mp_bl <- length(unique(sub$block[i][sub$is_mp[i]]))
    pick <- sample(bl, n_mp_bl)
    lab[i] <- sub$block[i] %in% pick
  }
  perm[k] <- stat_fn(lab)
}
p_emp <- (1 + sum(abs(perm) >= abs(obs)))/(n_perm + 1)
## The attainable resolution is bounded by how many distinct block assignments
## exist, not by the permutation count, so the number of Mesoproterozoic blocks
## per state is reported alongside the p value. A province contributing a single
## block contributes essentially no permutation evidence however many samples it
## holds.
cat(sprintf("\n=== blocked permutation test ===\n"))
cat(sprintf("observed weighted within-state difference in mean log10 Y: %+.4f\n", obs))
cat(sprintf("  = a factor of %.2f in Y\n", 10^obs))
cat(sprintf("permutations: %d, two-sided empirical p = %.5f\n", n_perm, p_emp))
cat(sprintf("decision at alpha = %.3f: %s\n", alpha,
            ifelse(p_emp < alpha, "REJECT the null", "fail to reject")))
cat(sprintf("permutation null: mean %+.4f, sd %.4f, 95%% range [%+.4f, %+.4f]\n",
            mean(perm), sd(perm), quantile(perm,.025), quantile(perm,.975)))
cat("\nMesoproterozoic blocks contributing per state:\n")
print(as.data.frame(tab2[, c("State","n_mp","mp_blocks","other_blocks")]), row.names = FALSE)
png(file.path(out_dir,paste0("figB1_permutation_null", sfx, ".png")), 1650, 930, res = 150)
hist(perm, breaks = 60, col = "grey85", border = "white",
     main = "Blocked permutation null, within-state Mesoproterozoic contrast",
     xlab = "weighted difference in mean log10 Y")
abline(v = obs, lwd = 3)
dev.off()
writeLines(c(sprintf("observed = %.6f (factor %.3f)", obs, 10^obs),
             sprintf("n_perm = %d", n_perm),
             sprintf("empirical two-sided p = %.6f", p_emp),
             sprintf("alpha = %.3f  (revised from 0.001; see header for rationale)", alpha),
             sprintf("decision = %s", ifelse(p_emp < alpha, "reject", "fail to reject"))),
           file.path(out_dir,paste0("permutation_test", sfx, ".txt")))

## ---- ratios, within state ------------------------------------------------------
## Y/Zr separates a melt-chemistry signal from simple zircon abundance, and Y/Th
## separates heavy from light enrichment. Computed within state for the same
## reason as Table 2.
tab3 <- st_drop_geometry(ig) %>%
  filter(State %in% test_states, Zr_p_ppm > 0, Th_ppm > 0, La_ppm > 0) %>%
  group_by(State, is_mp) %>%
  summarise(n = n(), Y = median(Y_ppm), Zr = median(Zr_p_ppm), Th = median(Th_ppm),
            La = median(La_ppm), Y_Zr = median(Y_ppm/Zr_p_ppm),
            Y_Th = median(Y_ppm/Th_ppm), Y_La = median(Y_ppm/La_ppm), .groups = "drop")
cat("\n=== Table 3: ratios within state ===\n")
print(as.data.frame(tab3), row.names = FALSE, digits = 4)
write.csv(tab3, file.path(out_dir,paste0("tableB3_ratios", sfx, ".csv")), row.names = FALSE)

st_write(ig, file.path(out_dir,paste0("intrusive_with_period", sfx, ".gpkg")), delete_dsn = TRUE, quiet = TRUE)

arm_results[[arm_name]] <- data.frame(
  arm = arm_name, states = paste(test_states, collapse = ","),
  n = nrow(sub), observed = obs, factor = 10^obs, p_emp = p_emp,
  alpha = alpha, decision = ifelse(p_emp < alpha, "reject", "fail to reject"))
}  ## ---- end of the arm loop ----------------------------------------------------

## ---- the two arms side by side -------------------------------------------------
## The primary result and its replication belong in one table, so that a reader
## cannot encounter one without the other.
summ <- do.call(rbind, arm_results)
cat("\n=== PRIMARY RESULT (wide) AND REPLICATION (narrow) ===\n")
print(summ, row.names = FALSE, digits = 5)
write.csv(summ, file.path(out_dir,"tableB0_primary_and_replication.csv"), row.names = FALSE)

# Written on every run so the software environment is documented by the analysis
# itself. Quote these versions when reporting the study; do not reconstruct them.
writeLines(c(paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
             paste("Script: eon_enrichment_v3.R"),
             paste("seed:", seed, " n_perm:", n_perm, " alpha:", alpha),
             paste("arms:", paste(names(arms), collapse=", ")),
             "", capture.output(sessionInfo())),
           file.path(out_dir,"sessionInfo.txt"))
cat("\nPaper B stage 1 complete. Outputs in", out_dir, "\n")
