% Abstract
\begin{abstract}
  A quick introduction to my favorite text editor \texttt{emacs}. Is
  it really only a text editor?
\end{abstract}

% Enable the following when using a report document style

%\tableofcontents
%\listoffigures
%\listoftables

\section{Introduction}
This document is mainly (and by now completely) based on the tutorials
in youtube by Baris Yuksel (\textit{b yuksel}),
\cite{emacs_video_01:URL, emacs_video_02:URL}, and the
\textit{Emacs-tutorial-for-Beginners} \cite{emacs_github:URL}, also by
Baris Yuksel.

In the following sections we present a set of commands that will make
you start talking \texttt{emacs} language. An \texttt{emacs} command
looks something like \texttt{C-x C-c}, in this example, the capital
\texttt{C} indicates to press and hold the \texttt{Control} key, the
lower case letters indicate to press the \texttt{x} and the \texttt{c}
keys on the keyboard. To introduce that command into \texttt{emacs}
you need to press and hold the \texttt{Control} key, then press the
\texttt{x} key, then press and hold the \texttt{Control} key and press
the \texttt{c} key, oh my gosh!!. You can also execute the command by
pressing and holding the \texttt{Control} key, then pressing the
\texttt{x} key followed by pressing the \texttt{c} key. Note that in
this case we pressed and hold the \texttt{Control} key and did not
released until once all the keys have been pressed.

You will find other commands having a capital \texttt{M}, this
indicates to press a \textit{Meta} key, which in most of the cases is
the \texttt{Alt} key. You can also use the \texttt{Esc} key, but it
may be kind of unconfortable while typing. Google for the meta key
corresponding to your system (Mac users).

\subsection{Is this your first time with \texttt{emacs}?}
I am not going to lie you, at the beginning you will find difficult
and weird to use the \texttt{emacs} commands, but once you get used to
them you will find that your typing speed has increased, and more
noticeable, you will no longer be using the mouse for typing your
documents. Well, in the first instance, why should we be using a mouse
for \textit{typing}?

Thus, do not panic if you find that you cannot copy and paste text,
you first need to learn the \texttt{emacs} language first so it
understand what you want.

Congratulations and welcome to the \texttt{emacs} world.

\section{Opening and quiting \texttt{emacs}}
\begin{itemize}
\item \texttt{emacs}: Opens \texttt{emacs} from a terminal.
\item \texttt{emacs -nw}: Opens \texttt{emacs} in a terminal. Use
  this if you have no access to X11 tools in your system.
\item \texttt{C-x C-c}: If the file has not been saved, it asks for
  saving it and then quits \texttt{emacs}.
\end{itemize}

\section{Visiting and saving files}
In \texttt{emacs}, opening or creating a new file is referred as
\textit{visiting} the file.
\begin{itemize}
\item \texttt{C-x C-f}: Opens a file, asks for the file name. If the
  filename does not exist then it creates a new file with the given
  name.
\item \texttt{C-x C-s}: Saves the file without a prompt.
\item \texttt{C-x s}: Saves all files with a prompt.
\item \texttt{C-s C-w}: Saves the file with a different name. Asks
  you for the name.
\end{itemize}

\subsection{Recovering files}
Let say you create a file called \texttt{my\_file.txt}, then emacs
automatically creates a file called \texttt{my\_file.txt\~} in the
same folder. This tilde(\texttt{~}) file is the previous version of
the file. Also, \texttt{emacs} has auto-save enabled by default, this
auto-save file is located in the same directory with the name
\texttt{\#my\_file.txt}. If for any reason (light cut) you quit
\texttt{emacs} without saving your file, you can recover it by opening
emacs and type \texttt{M-x recover-file}.

\section{Deleting text}
You can always use the \texttt{BackSpace} and \texttt{Delete} keys if
you want, but try thes ones and see if you still prefer those two old
fashioned keys.
\begin{itemize}
\item \texttt{C-d}: Deletes the letter at the cursor, same as
  \texttt{Delete}.
\item \texttt{M-d}: Deletes the word in front of the cursor, yes, the
  word. Think of this as a hungry \texttt{Delete}.
\item \texttt{C-k}: Deletes (and stores into the clipboard) the line
  in front of the cursor, yes, you read right, the complete
  line. Think of this as a very hungry \texttt{Delete}.
\end{itemize}

\section{Kill and Yank in \texttt{emacs}, or Cut/Copy/Paste}
In \texttt{emacs}, cutting a region or section of text is called
\textit{Kill} text, possibly because the text disappears from the
screen and is temporarly stored in the clipboard. \texttt{Emacs} use a
kind of \textit{stack} for the clipboard, you get access to the
history of elements stored in the clipboard by using a special command
for pasting (\texttt{M-y}).

\begin{itemize}
\item \texttt{C-space}: Starts marking/highligting a region. You can
  a large number of commands to the marked region, not just kill and
  yank.
\item \texttt{C-w}: (cut) Cuts this region into the clipboard
  (deletes the region and copies it to clipboard).
\item \texttt{C-k}: (cut) Kills/deletes the whole line, puts it into
  the clipboard.
\item \texttt{M-w}: (copy) Copies the selected region into the
  clipboard. Saving a region involves hitting C-space to start
  selecting, and then hitting \texttt{M-w} or \texttt{C-w} to copy or
  cut it into the clipboard, and then hitting \texttt{C-y} to paste
  it.
\item \texttt{C-y}: (paste) Pastes whatever is in the clipboard at
  the cursor. Subsequent \texttt{C-y}'s will keep on pasting.
\item \texttt{M-y}: (paste) Pastes whatever is in the clipboard at
  the cursor. Subsequent \texttt{M-y}'s will loop over the clipboard
  history.
\item \texttt{C-g}: Quits/cancels your command. If you dont like the
  region you are selecting, hit \texttt{C-g}.
\end{itemize}

\section{Moving around or cursor commands}
You can move around a buffer with need of the mouse, learn this few
commands and say bye bye to your mouse.
\begin{itemize}
\item \texttt{C-a}: Beginning of line.
\item \texttt{C-e}: End of line
\item \texttt{M->}: End of buffer
\item \texttt{M-<}: Beginning of buffer
\end{itemize}

\section{Thank gosh for undo and redo}
\texttt{Emacas} allows you to undo and redo, these are the commands:
\begin{itemize}
\item \texttt{C-/}: Undo
\item \texttt{C-g C-/}: Redo
\end{itemize}

\section{Frames, windows and buffers}
In \texttt{emacs}, each opened file gets a \textit{buffer}, there you
can edit the file and once you are done, save and close it. The buffer
used to visualise your file is then closed and released.

If you opened \texttt{emacs} in a terminal, you may want to visualise
more that one file at once, you can use \textit{windows} to split your
workspace and visualise more than one buffer at a time.

What you currently know as windows (not \texttt{emacs} windows), in
\texttt{emacs} languages are called \textit{frames}. Thus you can also
create new frames to visualise you files, note that inside a frame you
can have windows, and in each window a buffer.

\begin{itemize}
\item \texttt{C-x b}: Switches buffers, asks you which buffer to
  switch to.
\item \texttt{C-x C-b}: Switches buffers, but shows you the list of
  buffers in a new window.
\item \texttt{C-x o}: Move to other window. If you hit \texttt{C-x o}
  after showing the list of buffers then you can move over the name of
  the buffers, when you hit Enter over a buffer name then that buffer
  is opened in the current window.
\item \texttt{C-x 0}: Close the current window.
\item \texttt{C-x 1}: Close all other windows, and leave only the
  current one.
\item \texttt{C-x 2}: Make an horizontal cut to the current window to
  show a secondary window. 
\item \texttt{C-x 3}: Make a vertical cut and show a secondary
  window.
\end{itemize}

\section{Search mode}
One of functions that I use the most is the search function, in
\texttt{emacs} we can do different types of searchs/replacemetns in a
buffer.

\begin{itemize}
\item \texttt{C-s}: Searches \textit{forward} as you type. Jumps to
  the first instance after the cursor that matches what you
  typed. Other instances in the buffer are highlighted. You can go
  from one instance to the next one by hitting \texttt{C-s}
  again. When reaching the end of the buffer, then the search
  wraps-around starting from the beginning of the buffer. If you press
  \texttt{C-g} you can quit the search and return where you were.
\item \texttt{C-r}: Searches \textit{backwards} as you type. Jumps to
  the first instance before the cursor that matches what you
  typed. Other instances in the buffer are highlighted. You can go
  from one instance to the previous one by hitting \texttt{C-r}
  again. When reaching the beginning of the buffer, then the search
  wraps-around starting from the end of the buffer. If you press
  \texttt{C-g} you can quit the search and return where you were.
\item \texttt{M-\%}: Searches and replaces. You are asked to replace
  each term matching the search. If you hit \texttt{!} then it acst as
  replace all.
\item \texttt{M-C-\%}: Searches and replaces for regular
  expressions. You are asked to replace each term matching the
  search. If you hit \texttt{!} then it acst as replace all.
\item \texttt{M-s o}: Searches for regular expression and shows the
  matches in an \texttt{Occur} buffer, we can move into that buffer
  and hit Enter over the match to go to the specific place of the
  match.
\item \texttt{M-x grep}: greps a pattern in the files you specify and
  shows the results in a \texttt{grep} buffer. You can click or enter
  in the results to jump to the file and line where is the
  occurrance. The grep command uses this syntax: \texttt{grep -nH -e
    "string\_to\_search\_for" folder}, where the \texttt{-n} states to
  show lines numbers, the \texttt{H} is to show the file name and the
  \texttt{-e} is to indicate a pattern to search for. You can add
  \texttt{--colour} option to show colored output.
\item \texttt{M-x rgrep}: Recursive \texttt{grep} which searches in
  all the files and all the sub-directories in the given directory.
\end{itemize}

