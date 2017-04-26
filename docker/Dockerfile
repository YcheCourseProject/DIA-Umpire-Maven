FROM biocontainers/dia-umpire:1.4256

MAINTAINER Yulin CHE yche@cse.ust.hk

USER root

# install oh-my-zsh
RUN apt-get update && apt-get install -y \
    zsh \
    git-core && rm -rf /var/lib/apt/lists/*
RUN git clone https://github.com/robbyrussell/oh-my-zsh.git ~/.oh-my-zsh \
    && cp ~/.oh-my-zsh/templates/zshrc.zsh-template ~/.zshrc \
    && chsh -s /bin/zsh
