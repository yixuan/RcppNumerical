\name{NEWS}
\title{News for Package "RcppNumerical"}

\section{Changes in RcppNumerical version 0.4-0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Added one-dimensional numerical integration on infinite intervals,
            contributed by \href{https://github.com/rstub}{@rstub}
            (\href{https://github.com/yixuan/RcppNumerical/pull/11}{#11}).
    }
  }
}

\section{Changes in RcppNumerical version 0.3-3}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Updated LBFGS++ to the newest version.
      \item Fixed the \code{VignetteBuilder} entry in the \code{DESCRIPTION} file
            according to CRAN's policy.
    }
  }
}

\section{Changes in RcppNumerical version 0.3-2}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed a convergence test bug in LBFGS++.
      \item Registered native routines per CRAN's policy.
    }
  }
}

\section{Changes in RcppNumerical version 0.3-1}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed a bug of name hiding in C++.
    }
  }
}

\section{Changes in RcppNumerical version 0.3-0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Now using the LBFGS++ library for L-BFGS optimization.
      \item Added a vignette.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed an index-out-of-bound bug in one dimensional integration (#1),
            thanks to Greg Finak.
      \item Added virtual destructors for classes \code{Func}, \code{MFunc}, and
            \code{MFuncGrad}.
    }
  }
}

\section{Changes in RcppNumerical version 0.2-0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Added the Cuba library for multi-dimensional integration.
    }
  }
}

\section{Changes in RcppNumerical version 0.1-0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Initial release.
      \item Provides C++ interface for one-dimensional numerical integration
            and L-BFGS optimization.
      \item R function \code{fastLR()} for fast logistic regression using
            L-BFGS.
    }
  }
}
