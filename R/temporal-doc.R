#' Temporal assignment of uncertain recording dates
#'
#' Period-based assessment functions accept a `dateHandling` list for records
#' whose observation date is bounded by earliest and latest possible dates.
#'
#' @name temporal-assignment
#' @section Interface:
#' With `dateHandling = NULL`, the column named by `year` supplies exact calendar
#' years. Interval mode uses a named list such as:
#'
#' ```r
#' list(
#'   start = "startDate",
#'   end = "endDate",
#'   ambiguous = "maxOverlap",
#'   endOnly = "drop",
#'   startOnly = "drop",
#'   ties = "drop"
#' )
#' ```
#'
#' `start` and `end` name the earliest and latest possible observation-date
#' columns. Values may be numeric calendar years, `Date` values, four-digit year
#' strings, or ISO `YYYY-MM-DD` strings. Malformed values, impossible dates, and
#' starts after ends are errors. In interval mode, `year = NULL` is permitted.
#' If both bounds are missing, a usable exact year is the fallback. If either
#' bound is present, the interval information takes precedence over `year`.
#'
#' @section Assignment rules:
#' A bounded interval is assigned automatically when every calendar year it
#' touches is contained in one supplied period. Gaps in a period therefore
#' prevent automatic containment. Every retained record is assigned to at most
#' one period.
#'
#' For other bounded intervals, `ambiguous` is one of:
#'
#' - `"drop"` (default): exclude the record.
#' - `"start"`: use the period containing the earliest possible date.
#' - `"end"`: use the period containing the latest possible date.
#' - `"maxOverlap"`: use the period with greatest overlap. Year-only bounds use
#'   overlapping calendar-year counts; full dates use inclusive days within the
#'   actual years listed in each period.
#'
#' Maximum-overlap ties use `ties = "drop"` (default), `"earliest"`, or
#' `"latest"`. Earliest/latest choices are chronological and do not depend on
#' list order. A bounded interval with no overlap is excluded.
#'
#' An end-only record means "on or before" its bound. `endOnly = "drop"` is the
#' default; `"earliest"` assigns the chronologically earliest compatible period.
#' A start-only record means "on or after" its bound. `startOnly = "drop"` is the
#' default; `"latest"` assigns the chronologically latest compatible period.
#' These allocations are analyst assumptions, not information contained in the
#' record.
#'
#' @section Interpretation and diagnostics:
#' Interval allocation does not estimate the true date, duplicate records,
#' weight them fractionally, or imply continuous observation throughout a date
#' range. `maxOverlap` is a pragmatic allocation rule, not a probability model.
#' Preferential early or late allocation can change apparent temporal trends.
#'
#' Results retain their established named components. A
#' `temporalDiagnostics` attribute reports unambiguous and rule-based
#' assignments, exclusions by reason, unresolved ties, open-ended decisions,
#' outside-window records, missing temporal information, and final counts per
#' supplied period. Material interval exclusions or assumption-based allocations
#' produce one summary warning rather than one warning per record.
#'
#' Darwin Core `eventDate` can contain ISO dates or intervals. This interface
#' does not parse the complete Darwin Core interval syntax directly: split
#' bounded values into start and end columns first, retaining missing bounds for
#' open-ended records.
#'
#' @keywords internal
NULL
