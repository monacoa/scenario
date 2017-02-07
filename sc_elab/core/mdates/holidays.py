
import collections
import datetime

import dateutils

class OrderMapper:

    def __init__(self, items=None):
        self._map = dict()
        self._counter = 0
        if items:
            self.add(items)

    def add(self, items):
        for k in items:
            self._map[k] = self._counter
            self._counter += 1

    def __getitem__(self, item):
        return self._map[item]

    def __contains__(self, item):
        return item in self._map

# Names and abbreviated names of days of the week
# in English (not in locale, as definition of the calendars
# are written in English
day_names = [ 'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday' ]

MONDAY = 0
TUESDAY = 1
WEDNESDAY = 2
THURSDAY = 3
FRIDAY = 4
SATURDAY = 5
SUNDAY = 6

abbreviated_day_names = [ 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun' ]
month_names = ['', 'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August',
               'September', 'October', 'November', 'December']
month_name_order = OrderMapper(map(str.lower, month_names))
weekday_name_order = OrderMapper(map(str.lower, day_names))

def easter(year):
    """
    Calculate the date of Easter in the given year
    Adapted from Wikipedia: http://en.wikipedia.org/wiki/Computus
    """
    a = year % 19
    b = year // 100
    c = year % 100
    d = b // 4
    e = b % 4
    f = (b+8) // 25
    g = (b - f + 1) // 3
    h = (19*a+b-d-g+15) % 30
    i = c // 4
    k = c % 4
    L = (32+2*e+2*i-h-k) % 7
    m = (a + 11*h+22*L) // 451
    month = (h+L-7*m+114) // 31
    day = ((h+L-7*m+114) % 31)+1
    return datetime.datetime(year, month, day)

def holy_thursday(year):
    """ 
    For a given year determine a date of a Holy Thursday (Thursday before Easter)
    """
    return easter(year) - datetime.timedelta(days=3)

def good_friday(year):
    """ 
    For a given year determine a date of a Good Friday (Friday before Easter)
    """    
    return easter(year) - datetime.timedelta(days=2)

def easter_monday(year):
    """ 
    For a given year determine a date of an Easter Monday (Monday after Easter)
    """        
    return easter(year) + datetime.timedelta(days=1)

def ascension_thursday(year):
    """ 
    For a given year determine a date of an Ascension Thursday (Thursday 40 days after Easter)
    """            
    return easter(year) + datetime.timedelta(days=39)

def pentecost(year):
    """ 
    For a given year determine a date of Pentecost (Sunday 7 weeks from Easter)
    """            
    return easter(year) + datetime.timedelta(days=49)

def whit_monday(year):
    """ 
    For a given year determine a date of a Whit Monday (50 days from Easter)
    """                
    return easter(year) + datetime.timedelta(days=50)

def trinity_sunday(year):
    """ 
    For a given year determine a date of a Trinity (Sunday 8 weeks from Easter)
    """                    
    return easter(year) + datetime.timedelta(days=56)

def corpus_christi_thursday(year):
    """ 
    For a given year determine a date of a Trinity (Sunday 8 weeks from Easter)
    """                        
    return easter(year) + datetime.timedelta(days=60)

def weekday_on_or_before(dt, weekday):
    """ Find weekday happening on or before a given date
    
    This method is used to determine the date of some holidays, e.g.
    Victoria day in Canada is celebrated on the Monday on or 
    preceding 24th of May
    
    Parameters
    ----------
    dt : datetime.date or datetime.datetime
    weekday: weekday number (0 - Monday, 6 - Sunday)
    
    Returns
    -------
    datetime.date or datetime.datetime with a date of a holiday
    """
    
    end_wd = dt.weekday()
    delta = end_wd - weekday
    if delta < 0:
        delta += 7
    return dt - datetime.timedelta(days=delta)
    

def weekday_on_or_after(dt, weekday):
    """ Find weekday happening on or after a given date
    
    This method is used to determine the date of some holidays, e.g.
    Victoria Day in Canada is celebrated on Monday on or before 24th of May
    
    Parameters
    ----------
    dt : datetime.date or datetime.datetime
    weekday: weekday number (0 - Monday, 6 - Sunday)
    
    Returns
    -------
    datetime.date or datetime.datetime with a date of a holiday
    """
    
    start_wd = dt.weekday()
    delta = weekday - start_wd
    if delta < 0:
        delta += 7
    return dt + datetime.timedelta(days=delta)
    
def midsummer_eve(year):
    """ Celebrated in Finland and Sweden: Friday between June 18-24 inclusive
    """
    dt_start = datetime.datetime(year, 6, 18)
    return weekday_on_or_after(dt_start, FRIDAY)

def victoria_day(year):
    """ Victoria day in Canada - Monday on or preceding 24th of May
    """
    dt_end = datetime.datetime(year, 5, 24)
    return weekday_on_or_before(dt_end, MONDAY)

idiosyncratic_holidays = dict({
    'holy thursday': holy_thursday,
    'good friday': good_friday,
    'easter monday': easter_monday,
    'ascension thursday': ascension_thursday,
    'pentecost': pentecost,
    'whit monday': whit_monday,
    'trinity sunday': trinity_sunday,
    'corpus christi thursday': corpus_christi_thursday,
    'midsummer eve': midsummer_eve,
    'victoria day': victoria_day
    })

class Calendar:
    def __init__(self):
        self._weekdays = [True] * dateutils.DAYS_IN_WEEK
        self._dated_holidays = collections.defaultdict(dict)
        self._numbered_weekday_holidays = collections.defaultdict(set)
        self._holiday_cache = collections.defaultdict(set)
        self._idiosyncratic_years_cached = dict()
        self._idiosyncratic = dict()

    def add_holiday(self, name, date_description, **kwargs):
        """ Add holiday to the calendar

        Parameters
        ----------
        name: name of the holiday (use 'weekend' for weekend days as it has special
              meaning for holiday moving rules)
        date_description: date description in one of the following forms:
            * A name of the day of the week. Most probably you must specify
              'name' parameter as 'weekend' in this case
            * a specific date in "month day" string (e.g. "July 1st"
            * specific weekday in the month, e.g. "2nd Thursday in June" or
              "last Monday in May"
            * a name of holiday for which add_holiday supports a special calculation,
              e.g. 'Pentecost' or 'Victoria Day'
        Keyword parameters:
        move: can be 'next' or 'closest'. If holiday falls on a weekend (days of week
          specified with 'weekend' in the previous 'add_holiday()' call, then holiday is
          moved to the next available day that was not designated as weekend or holiday.
          'closest' first try to move holiday to both the available previous day and available
          following day and picks the one closest to the actual holiday date. In that case, if you
          have Saturdays and Sundays as weekend days, holidays happening on Saturday will be moved
          to Friday and those happening on Sunday will be moved to Monday.

        """
        if 'move' in kwargs:
            move = kwargs['move']
        else:
            move = None
        # go through days of the week
        if date_description.lower() in idiosyncratic_holidays.keys():
            self._idiosyncratic[date_description.lower()] = (name, move)
        for day_name_idx in range(len(day_names)):
            day_name = day_names[day_name_idx]
            if date_description.lower()==day_name.lower():
                self._weekdays[day_name_idx] = False
        desc_parts = date_description.split(' ')
        if len(desc_parts)==2:
            month_name = desc_parts[0].lower()
            if month_name in month_name_order:
                daystr= desc_parts[1].lower()
                if daystr.endswith('st') or daystr.endswith('nd') or daystr.endswith('rd') or daystr.endswith('th'):
                    daystr = daystr[:-2]
                day_num = int(daystr)
                month_num = month_name_order[month_name]
                self._dated_holidays[month_num][day_num] = (name, move)
        if len(desc_parts)==4:
            # nth weekday in month
            month_num = month_name_order[desc_parts[3].lower()]
            weekday_num = weekday_name_order[desc_parts[1].lower()]
            order_str = desc_parts[0]
            if order_str.lower() == 'last':
                order = -1
            else:
                if order_str.endswith('st') or order_str.endswith('nd') or order_str.endswith('rd') or order_str.endswith('th'):
                    order_str = order_str[:-2]
                order = int(order_str)
            self._numbered_weekday_holidays[month_num].add( (order, weekday_num, name) )

    def _is_idiosyncratic(self, dt):
        """ check if date is one of idiosyncratic holidays
        """
        Y = dt.year
        if not Y in self._idiosyncratic_years_cached:
            # compute idiosyncratic holidays for year Y
            id_table = dict()
            for id in self._idiosyncratic.keys():
                name, move = self._idiosyncratic[id]
                id_table[idiosyncratic_holidays[id](Y)] = (name, move)
            self._idiosyncratic_years_cached[Y] = id_table
        if dt in self._idiosyncratic_years_cached[Y]:
            name, move = self._idiosyncratic_years_cached[Y][dt]
            return True, name, move
        else:
            return False, None, None

    def is_holiday(self, dt):
        """ Check if specific date is holiday
        """
        dt = dateutils.asdatetime(dt)
        year = dt.year
        if not year in self._holiday_cache:
            # create year cache
            self._holiday_cache[year] = set()
            current = datetime.datetime(year, 1, 1)
            end_year = datetime.datetime(year, 12, 31)
            while current <= end_year:
                b_holiday, name, move_day = self._verify_holiday(current)
                if b_holiday:
                    self._holiday_cache[year].add(current)
                if move_day is not None:
                    self._holiday_cache[year].add(move_day)
                current = current + datetime.timedelta(days = 1)
        return dt in self._holiday_cache[year]

    def _move_holiday(self, dt, move):
        next_day = dt + datetime.timedelta(days = 1)
        while not self._weekdays[next_day.weekday()] or self._verify_holiday(next_day)[0]:
            next_day += datetime.timedelta(days = 1)
        if move == 'closest':
            prev_day = dt - datetime.timedelta(days = 1)
            while not self._weekdays[prev_day.weekday()] or self._verify_holiday(prev_day)[0]:
                prev_day -= datetime.timedelta(days = 1)
            delta_prev = dt - prev_day
            delta_next = next_day - dt
            if delta_prev < delta_next:
                result = prev_day
            else:
                result = next_day
        else:
            # add next day to the list of moves
            result = next_day
        return result


    def _verify_holiday(self, dt):
        move_day = None
        weekday = dt.weekday()

        is_weekend = not self._weekdays[weekday]

        idiosyncratic, name, move = self._is_idiosyncratic(dt)
        # No move rules for idiosyncratic holidays?
        if idiosyncratic:
            return True, name, None

        if dt.month in self._dated_holidays:
            holiday_day_descriptions = self._dated_holidays[dt.month]
            if dt.day in holiday_day_descriptions:
                name, move = holiday_day_descriptions[dt.day]
                if move is not None and is_weekend:
                    move_day = self._move_holiday(dt, move)
                return True, name, move_day

        if len(self._numbered_weekday_holidays[dt.month]) > 0:
            for tup in self._numbered_weekday_holidays[dt.month]:
                order, weekday, name = tup
                # n-th day of the month holidays are usually not moved
                # as they always happen on a particular day of the week
                if dt.weekday() == weekday:
                    if order > 0 and (((dt.day-1)//dateutils.DAYS_IN_WEEK+1) == order):
                        return True, name, None
                    if order == -1:
                        # last weekday of the month
                        eom_day = dateutils.eom(dt.year, dt.month)
                        if dt == weekday_on_or_before(eom_day, weekday):
                            return True, name, None

        # if everything else in other categories does not work - check if it is weekend
        if is_weekend:
            return True, 'weekend', None

        return False, '', move_day

def get_calendar(calendar_code):
    cl = Calendar()
    calendar_code = calendar_code.lower()
    if calendar_code == 'us' or calendar_code == 'united states':
        cl.add_holiday("weekend", "Saturday")
        cl.add_holiday("weekend", "Sunday")
        cl.add_holiday("New Year's Day", "January 1st", move='closest')
        cl.add_holiday("Martin Luther King's birthday", "3rd Monday in January")
        cl.add_holiday("Presidents' day", "3rd Monday in February")
        cl.add_holiday("Independence Day", "July 4th", move='closest')
        cl.add_holiday("Labor Day", "1st Monday in September")
        cl.add_holiday("Columbus Day", "2nd Monday in October")
        cl.add_holiday("Veterans Day", "November 11th", move='closest')
        cl.add_holiday("Thanksgiving", "4th Thursday in November")
        cl.add_holiday("Christmas", "December 25th", move='closest')
    elif calendar_code == "ca" or calendar_code == "canada":
        cl.add_holiday("weekend", "Saturday")
        cl.add_holiday("weekend", "Sunday")
        cl.add_holiday("New Year's Day", "January 1st", move='next')
        cl.add_holiday("Good Friday", "Good Friday")
        cl.add_holiday("Easter Monday", "Easter Monday")
        cl.add_holiday("Victoria Day", "Victoria Day")
        cl.add_holiday("Canada Day", "July 1st", move='next')
        cl.add_holiday("Civic Holiday", "1st Monday in August")
        cl.add_holiday("Labor Day", "1st Monday in September")
        cl.add_holiday("Thanksgiving", "2nd Monday in October")
        cl.add_holiday("Remembrance Day", "November 11th")
        cl.add_holiday("Christmas", "December 25th", move='next')
        cl.add_holiday("Boxing Day", "December 26th", move='next')
    elif calendar_code == 'de':
        cl.add_holiday("weekend", "Saturday")
        cl.add_holiday("weekend", "Sunday")
        cl.add_holiday("New Year's Day", "January 1st", move='next')
        cl.add_holiday("Good Friday", "Good Friday")
        cl.add_holiday("Easter Monday", "Easter Monday")
        cl.add_holiday("Ascension Thursday", "Ascension Thursday")
        cl.add_holiday("Whit Monday", "Whit Monday")
        cl.add_holiday("Corpus Christi Thursday", "Corpus Christi Thursday")
        cl.add_holiday("Labour Day", "May 1st")
        cl.add_holiday("National Day", "October 3rd")
        cl.add_holiday("Christmas Eve", "December 24th")
        cl.add_holiday("Christmas", "December 25th")
        cl.add_holiday("Boxing Day", "December 26th")
        cl.add_holiday("New Year's Eve", "December 31st")
    elif calendar_code == 'de.frankfurt' or calendar_code == 'de.xetra' or calendar_code == 'de.eurex':
        cl.add_holiday("weekend", "Saturday")
        cl.add_holiday("weekend", "Sunday")
        cl.add_holiday("New Year's Day", "January 1st", move='next')
        cl.add_holiday("Good Friday", "Good Friday")
        cl.add_holiday("Easter Monday", "Easter Monday")
        cl.add_holiday("Labour Day", "May 1st")
        cl.add_holiday("Christmas Eve", "December 24th")
        cl.add_holiday("Christmas", "December 25th")
        cl.add_holiday("Boxing Day", "December 26th")
        cl.add_holiday("New Year's Eve", "December 31st")
    elif calendar_code == 'uk':
        cl.add_holiday("weekend", "Saturday")
        cl.add_holiday("weekend", "Sunday")
        cl.add_holiday("New Year's Day", "January 1st", move='next')
        cl.add_holiday("Good Friday", "Good Friday")
        cl.add_holiday("Easter Monday", "Easter Monday")
        cl.add_holiday("Early May Bank Holiday", "1st Monday in May")
        cl.add_holiday("Spring Bank Holiday", "last Monday in May")
        cl.add_holiday("Summer Bank Holiday", "last Monday in August")
        cl.add_holiday("Christmas", "December 25th", move='next')
        cl.add_holiday("Boxing Day", "December 26th", move='next')
    else:
        raise ValueError('unknown calendar code \'%s\'' % calendar_code)
    return cl
