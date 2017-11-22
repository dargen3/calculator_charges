from subprocess import check_call, DEVNULL


def alarm():
    try:
        check_call(["play", "modules/alarm/alarm.mp3"], stdout=DEVNULL, stderr=DEVNULL)
    except KeyboardInterrupt:
        pass
    except CalledProcessError:
        print(colored("You cannot run mp3 by bash command play!", "red"))
