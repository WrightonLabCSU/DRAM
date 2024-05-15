  [75] "blaze7B-rt-cf-qt_MH_bin-2_k121_1035912_15" 
  [76] "blaze11-rt-cf-qt_MH_bin-35_k121_167870_2"  
  [77] "blaze23-rt-cf-qt_MH_bin-8_k121_501691_19"  
  [78] "blaze7B-rt-cf-qt_MH_bin-29_k121_1242206_10"
  [79] "blaze22-rt-cf-qt_MH_bin-7_k121_809961_9"   

Command error:
      rotate
  
  Average angle change [1] 0.147309399409057
  Average angle change [2] 0.155623744117933
  Average angle change [3] 0.126625980724506
  Average angle change [4] 0.16382246464418
  Average angle change [5] 0.12569584526505
  Error in `geom_segment2()`:
  ! Problem while computing aesthetics.
  ℹ Error occurred in the 6th layer.
  Caused by error:
  ! object 'node' not found
  Backtrace:
       ▆
    1. ├─ggplot2::ggsave(output_pdf, plot = p, width = 20, height = 20)
    2. │ ├─grid::grid.draw(plot)
    3. │ └─ggplot2:::grid.draw.ggplot(plot)
    4. │   ├─base::print(x)
    5. │   └─ggplot2:::print.ggplot(x)
    6. │     ├─ggplot2::ggplot_build(x)
    7. │     └─ggplot2:::ggplot_build.ggplot(x)
    8. │       └─ggplot2:::by_layer(...)
    9. │         ├─rlang::try_fetch(...)
   10. │         │ ├─base::tryCatch(...)
   11. │         │ │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
   12. │         │ │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
   13. │         │ │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
   14. │         │ └─base::withCallingHandlers(...)
   15. │         └─ggplot2 (local) f(l = layers[[i]], d = data[[i]])
   16. │           └─l$compute_aesthetics(d, plot)
   17. │             └─ggplot2 (local) compute_aesthetics(..., self = self)
   18. │               └─base::lapply(aesthetics, eval_tidy, data = data, env = env)
   19. │                 └─rlang (local) FUN(X[[i]], ...)
   20. └─base::.handleSimpleError(`<fn>`, "object 'node' not found", base::quote(NULL))
   21.   └─rlang (local) h(simpleError(msg, call))
   22.     └─handlers[[1L]](cnd)
   23.       └─cli::cli_abort(...)
   24.         └─rlang::abort(...)
  Warning messages:
  1: In max(x, na.rm = TRUE) :
    no non-missing arguments to max; returning -Inf
  2: In min(x, na.rm = na.rm) :
    no non-missing arguments to min; returning Inf
  3: In max(x, na.rm = na.rm) :
    no non-missing arguments to max; returning -Inf
  4: In min(x, na.rm = na.rm) :
    no non-missing arguments to min; returning Inf
  5: In max(x, na.rm = na.rm) :
    no non-missing arguments to max; returning -Inf
  Execution halted
