from subprocess import check_call

def alarm():
    try:
        check_call("play modules/alarm/alarm.mp3 2> /dev/null", shell=True)
    except KeyboardInterrupt:
        pass
    except CalledProcessError:
        print(colored("You cannot run mp3 by bash command play!", "red"))