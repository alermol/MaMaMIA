# Karyotype tooltip reference

The karyotype tooltip is a call card for one introgression call. It answers four
questions in a fixed order:

1. What is the call?
2. Where was it placed?
3. How strong is the evidence?
4. How much sequence changed?

The card is shown when the pointer is over a chromosome bar in the
**Karyotype** tab. The content follows the current dashboard filters; a call
that is hidden from the call table does not appear in the tooltip. The
**Tooltip help** button beside the Heatmap/Calls switch opens an in-app
field guide covering the same fields, interpretation patterns and caveats.

## Card anatomy

A call card has up to six blocks:

```text
header
placement
evidence
dosage
geometry
position (heatmap mode only)
footer
```

The geometry block is omitted when no geometry field is available. Every other
block is present for a call card. Fields that are missing or `NA` are omitted
inside a block, except the dosage `Copies` row, where `unavailable` is useful
information.

## Header

| Element | Source | Meaning |
|---|---|---|
| Class title | `class` | `substitution`, `addition`, `deletion`, `complex` or `other`. The class colour is also the card accent. |
| Status pill | `class_status`; `class_confidence` is mapped to the same levels when the status is absent | `confident`, `ambiguous` or `needs_review`. |
| Context line | `sample`, `group`, caller coordinates | Sample name, chromosome group and donor span. |

The status pill describes confidence in the class, not confidence in the
recipient placement. A call can be `confident` as a substitution but still have
an `ambiguous` assignment.

## Placement

| Row | Source fields | Meaning |
|---|---|---|
| `Recipient` | `recipient_chr_name`, `recipient_group`, `recipient_start_bp`, `recipient_end_bp` | Assigned recipient contig and interval. Shows `no host` when `recipient_assignment` is `no_recipient_loss`. |
| `Donor` | `donor_chr_name`, `donor_group`, `donor_start_bp`, `donor_end_bp` | Donor contig and interval. |
| `Relation` | `cross_group`, `recipient_assignment` | `same group`, `cross-group`, or `not placed`. |
| `Assignment` | `recipient_assignment`, `recipient_margin` | Placement state and, when finite, the margin over the next candidate. |
| `Group class` | `class_group`, `class` | Shown only when assignment changed the class. For example, `addition → substitution`. |
| `Alternative` | `recipient_alternative` | Shown only when the assignment is `ambiguous` and another host exists. |

Assignment labels:

| Value | Label |
|---|---|
| `assigned` | `assigned` |
| `cross_group` | `cross-group candidate` |
| `ambiguous` | `ambiguous` |
| `no_recipient_loss` | `no host` |
| missing | `group only` |

`cross-group` means the best host is outside the call's chromosome group. It is
a translocation candidate, not a confirmed placement. `no host` means the call
has no recipient deletion to anchor it.

## Evidence

The evidence block starts with four metric cards:

| Metric | Source | Meaning |
|---|---|---|
| `Score` | `score_mean` | Paired donor/recipient score for the segment. |
| `Donor` | `donor_mean` | Mean donor-side residual over the segment. |
| `Recipient` | `recipient_mean` | Mean recipient-side residual over the segment. |
| `Total` | `total_mean` | Relative total-depth residual. |

Values are shown with two decimals. Positive values are green, negative values
red, and missing values muted. A positive recipient residual means recipient
gain; a negative recipient residual means recipient loss.

The metric cards are followed by zero or more chips. A chip appears only when
its condition holds.

| Chip | Source | Condition |
|---|---|---|
| `class high`, `class moderate`, `class low` | `class_confidence` | Always when the column is present. |
| `donor`, `recipient`, `donor+recipient`, `none` | `class_evidence` | Always when the column is present. Compound calls append the reason, for example `donor+recipient; total imbalance`. |
| `balanced`, `gain`, `loss`, `compound`, `neutral` | `copy_state` | Always when the column is present. |
| `total gain`, `total loss` | `total_evidence` | Only when total evidence is not `neutral`. |
| `components heterogeneous` | `component_status` | Only when the value is `heterogeneous`. |
| `low mappability` | `mappability_status` | Only when the value is `low`. |
| `not passing` | `pass` | Only when the call is not a passing call. |

`copy_state` is a hypothesis about the copy configuration, not a measured
absolute copy number. `components heterogeneous` means the component spread
exceeds the 2-log cutoff or more than one donor component is supported, so the
group-level score may hide a component-specific event. `low mappability` means
the BAM-derived high-MAPQ fraction is below 0.5.

## Dosage

| Row | Source | Meaning |
|---|---|---|
| `Copies` | `affected_copies`, `recipient_copy_number` | Affected copies, for example `2 of 2 copies`. Shows `unavailable` when no copy count is estimated. |
| `Estimated fraction` | `affected_fraction` | Affected fraction of recipient copies, as a percentage. |
| `Estimated donor copies` | `donor_copies` | Relative donor-copy estimate when finite. |
| `Confidence` | `dosage_confidence` | `high`, `moderate`, `low` or `unavailable`. |

Dosage is relative to the registry copy number. It is not an absolute
copy-number measurement. The fraction and donor-copy rows are labelled
`Estimated` because they are continuous estimates; only `Copies` is rounded to
an integer. For an addition, `Copies` comes from the donor-copy estimate. For a
deletion without donor support, it comes from the recipient loss fraction at
lower confidence.

## Geometry

| Row | Source | Meaning |
|---|---|---|
| `Donor width` | `width_bp` | Width of the donor segment. |
| `Recipient run` | `recipient_deletion_width_bp` | Width of the recipient deletion run used for assignment. |
| `Overlap` | `recipient_width_overlap` | Reciprocal overlap of donor and recipient intervals in contig-fraction space. |
| `Width ratio` | `recipient_fraction_ratio` | Donor fraction divided by recipient-run fraction. |

The whole block is omitted when none of these fields is finite. Width ratio and
overlap are the values used to decide whether a donor/recipient pair has
compatible geometry. A width ratio far from 1, or an overlap below 0.5, can turn
a substitution into `complex`.

## Position block

The position block appears only in **Heatmap** mode. It shows:

- `Position`: hovered position along the chromosome, formatted as bp, kb or Mb.
- `Score`: the heatmap score at that relative position.

The block is added after the call card. In **Calls** mode the card ends with the
footer instead.

## Footer

The footer always says:

```text
Click to filter this group
```

Clicking the chromosome sets the chromosome-group filter. Clicking the
karyotype background clears it.

## No-call card

When the pointer is not over a call, the tooltip becomes a chromosome card.

```text
RECIPIENT CHROMOSOME             [ no call ]
example-sample · group 1 · 1.9 kb

CHROMOSOME
Component  A
State      recipient chromosome
Homolog    1 of 2
Length     20.0 kb

CALL
No assigned donor segment at this position.
```

The title depends on the chromosome state:

| State | Header title | `State` row |
|---|---|---|
| Normal recipient | `recipient chromosome` | `recipient chromosome` |
| Deleted recipient | `deleted recipient` | `deleted recipient (no coverage)` |
| Added donor | `added donor` | `added donor chromosome` |

The `Homolog` row appears only when the chromosome has more than one configured
copy. The `CALL` section appears only in Calls mode. In Heatmap mode the
Position block is shown instead.

## Formatting rules

- Genomic distances use `Mb` at 1,000,000 bp or more, `kb` at 1,000 bp or more,
  and `bp` below that.
- Signed residuals and margins use two decimals and an explicit `+` for
  positive values.
- Percentages are rounded to whole numbers.
- Missing numeric values show `NA` in metric cards and are omitted from rows.
- The dosage `Copies` row is the exception: it shows `unavailable` when no copy
  count exists.
- The card is 340 px wide, capped to the viewport, and its height is capped to
  the viewport.
- The tooltip is fixed-positioned, flipped left or up near the viewport edge,
  and clamped inside the karyotype area.
- `pointer-events: none` keeps the tooltip from blocking chromosome hover.

## Visual encoding

- The left border and class title use the class colour from the karyotype
  palette.
- The status pill uses green for `confident`, amber for `ambiguous`, red for
  `needs_review`, and grey for `no call`.
- Metric cards use a four-column grid.
- Chips wrap below the metrics.
- Section labels are uppercase and separated by rules.
- Placement and geometry use a two-column label/value grid.

## Data flow

```text
BamAnalysisResult$calls
        │
        ▼
karyotype_payload()
  attaches calls to recipient or donor chromosomes
        │
        ▼
chromosome$calls
        │
        ▼
describe() in mamamiaKaryotypeCard.js
        │
        ▼
HTML card + mamamiaKaryotypeCard.css
```

`karyotype_payload()` in `R/karyotype.R` builds one payload entry per chromosome
copy. The payload contains one `calls` list per chromosome. Each call entry
carries the fields listed above. The JavaScript `describe()` function finds the
call covering the hovered fraction and renders the card.

Recipient bars include calls that have a recipient host. Calls with
`recipient_assignment = "no_recipient_loss"` are excluded from recipient bars
but can still appear on added donor chromosomes.

## Reading examples

| Card state | Interpretation |
|---|---|
| `substitution`, `balanced`, `assigned` | Donor sequence replaces recipient sequence with a neutral total depth. |
| `addition`, `gain`, `no host` or `assigned` | Donor sequence is gained without a matching recipient loss. |
| `deletion`, `loss`, `assigned` | Recipient sequence is lost. Donor support is absent or weak. |
| `complex`, `compound`, `needs_review` | Total depth or geometry is inconsistent with a simple substitution. |
| `ambiguous` with `Alternative` | A host exists, but the margin over the next candidate is small. |
| `cross-group candidate` | The best host is outside the call's group. Treat it as a translocation candidate. |
| `not passing` chip | The segment did not pass calibrated detection. It is visible only when the dashboard shows non-passing calls. |

## Limitations

- The card describes group-level evidence. A component-specific event appears
  only as `components heterogeneous`; it is not split into a separate card.
- Total depth is relative. The card cannot report an absolute total copy
  number.
- Dosage is relative to the registry copy number and depends on breadth.
- `low mappability` is a MAPQ-based proxy. It is not a full mappability track.
- Homolog identity is available only when the registry declares separate
  homolog contigs. Collapsed homologs are reported as a copy set.
- `class_status` and `recipient_assignment` answer different questions. Read
  both before treating a call as a confirmed introgression.

## Implementation and tests

| Concern | File |
|---|---|
| Payload construction | `R/karyotype.R` |
| Card renderer | `inst/htmlwidgets/mamamiaKaryotypeCard.js` |
| Card layout | `inst/htmlwidgets/mamamiaKaryotypeCard.css` |
| Payload contract test | `tests/testthat/test-karyotype.R` |

The widget assets use the name `mamamiaKaryotypeCard` and the CSS dependency
version `0.1.0`, so a browser cache cannot keep serving the previous flat
tooltip after an update. JavaScript syntax is checked with `node --check`.
