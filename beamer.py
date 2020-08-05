import os
import getpass

class BeamerPresentation(object):
    def __init__(self,name,title=None,subtitle=''):
        os.makedirs(name)
        self.f = open(os.path.join(name,name+'.tex'),'w')
        self.title = name if title is None else title
        self.subtitle = subtitle
        self.header()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.finish()

    def header(self):
        author = getpass.getuser()
        self.f.write(r"""%%% Generated with compare_image_dirs.py %%%
\documentclass{beamer}
\usetheme{default}

\title{Image Comparison}
\subtitle{"""+self.subtitle+r"""}
\author{"""+author+r"""}
\date{\today}
\begin{document}

\begin{frame}[plain]
    \titlepage
\end{frame}

""")

    def add_frame_2images(self,img1,img2,title='',subtitles=['','']):
        self.f.write(r'\begin{frame}{'+title+r"""}
\begin{columns}
    \begin{column}{0.5\textwidth}
        \begin{center}
        """+subtitles[0]+r"""
        \includegraphics[width=\textwidth]{"""+img1+r"""}
        \end{center}
    \end{column}
    \begin{column}{0.5\textwidth}
        \begin{center}
        """+subtitles[1]+r"""
        \includegraphics[width=\textwidth]{"""+img2+r"""}
        \end{center}
    \end{column}
\end{columns}
\end{frame}

""")

    def finish(self):
        self.f.write(r'\end{document}')
        self.f.close()

