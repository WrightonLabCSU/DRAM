# Rule Expression Grammar and Parser

DRAM uses a small, explicit expression language to define **rules**, **conditions**, and **logical combinations** of features.
This page documents the grammar, available operators and functions, and an important design decision regarding **boolean operator ambiguity**.

---

## Overview

Rule expressions are used throughout DRAM to:

* Summarize (distillate) annotated genes into a summarized sheet (see distill_metals.tsv for examples)
* Evaluate complex rules for gene traits (see trait_rules.tsv for examples)
* Build visualizations (see dram_viz for examples)

The grammar is intentionally **restrictive** to avoid ambiguous interpretations that can easily arise in boolean logic.

---

## Expression Components

The lowest-level building blocks of expressions are genes:
        
```text
KXXXXX (KEGG Orthology ID)
```

These can be evaluated indivudally as presence/absence in your annotation sheet and then combined using boolean logic on other genes to get totaly trait presence/absence or scores.

### Boolean Operators

| Operator | Meaning     | Notes                                         |
| -------- | ----------- | --------------------------------------------- |
| `&`      | Logical AND | May only be chained with other `&` operators  |
| `\|`     | Logical OR  | May only be chained with other `\|` operators |

❗ Mixing `&` and `\|` **requires explicit grouping**.


### Grouping

Square brackets are used to explicitly group expressions:

```text
[KXXXXX & KXXXXX] | KXXXXX
```

Brackets may contain **any valid expression**, including nested boolean logic.

### Step Expressions

Comma-separated expressions define **step rules**:

```text
KXXXXX & KXXXXX, KXXXXX | KXXXXX
```

These steps represent a sequence or a pathway of evaluations. This can be fed into functions that operate on multiple steps, such as looking for a percentage of steps satisfied.

---

## Core Design Principle: No Implicit Boolean Precedence

DRAM **does not allow mixing `|` (OR) and `&` (AND) operators without explicit grouping**.

This means expressions like:

```text
A | B | C & D
```

are **invalid** and will result in a parsing error.

Instead, users must write:

```text
A | B | [C & D]
```

or:

```text
[A | B | C] & D
```

### Why this matters

In many languages, `&` has higher precedence than `|`, but this is:

* Often misunderstood
* Easy to misread in complex biological rules
* A frequent source of subtle bugs

DRAM therefore **disallows implicit precedence** and requires explicit grouping with brackets (`[...]`) whenever boolean operators are mixed.

---

## Functions

Functions operate on expressions or values and may appear anywhere an expression is allowed.

### General Form

```text
function_name(arg1, arg2, ...)
```

Arguments may be:

* Expressions
* Numbers
* Strings
* Identifiers

---

### Functions

- not
- percent
- at_least
- column_count_values


#### not

Negates a boolean expression.

```text
not(A & B)
```

#### percent

Calculates percentage of steps satisfied and takes in a threshold value. Returns true if the percentage of satisfied steps meets or exceeds the threshold.

```text
percent(60, [A & B, C | D, E])
```

#### at_least

Checks if at least a certain number of steps are satisfied.

```text
at_least(2, [A & B, C | D, E])
```

#### column_count_values


```text
/**
 * Counts the number of values in a specified column and evaluates conditions based on the provided operations and thresholds.
 *
 * @param column_name The name of the column to be analyzed.
 * @param value_op The operation to be used for comparing the column values against the value_threshold (e.g., "gt", "ge", "lt", "le", "eq", "ne").
 * @param value_threshold The threshold value to compare the column values against using value_op.
 * @param count_op The operation to be used for comparing the counted values against the count_threshold (e.g., "gt", "ge", "lt", "le", "eq", "ne").
 * @param count_threshold The threshold value to compare the counted values against using count_op.
 * @return boolean Returns true if the conditions defined by value_op and count_op are satisfied; otherwise, returns false.
 */

column_count_values("column_name", "value_op", value_threshold, "count_op", count_threshold)
```

---

### Piped Function Usage

Filter functions may be used with the pipe operator:

```text
not(filter_contains(kegg_description,"nitrate reductase")) -> column_count_values(heme_regulatory_motif_count,ge,4,ge,3)
```

These filter functions can be used to prefilter data before applying counting or other operations. Filter functions are under the `filter_` namespace.

- filter_contains
- filter_compare

#### filter_contains

Checks if a specified column contains a given substring.

```text
filter_contains(column_name, substring)
```

#### filter_compare

Compares values in a specified column against a threshold using a given operation.

```text
filter_compare(column_name, operation, threshold)
```

---

## Valid vs Invalid Examples

### ✅ Valid

```text
A | B | C
A & B & C
[A & B] | C
A | [B & C]
filter_contains(col, "substring A") -> column_count_values(col,ge,2,ge,1)
```

---

### ❌ Invalid

```text
A | B | C & D
A & B | C
A | B & C | D
[A, B, C] -> percent(50)
```
